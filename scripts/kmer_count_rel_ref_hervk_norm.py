#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys
import pysam
import numpy as np
import log,traceback


cigar_op={'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
cigar_ref_retain={'M', 'D', 'N', '=', 'X'}
cigar_read_retain={'M', 'I', '=', 'X'}


def calc_depth_coeff(args, params):
    if not args.b is None:
        stdout=pysam.depth(args.b, '-a')
    else:
        stdout=pysam.depth(args.c, '-a', '--reference', args.fa)
    stdout=stdout.strip().split('\n')
    global coeff
    coeff={}
    d={}
    for line in stdout:
        ls=line.split()
        if ls[0] == args.refseq_id:
            d[int(ls[1]) - 1]=int(ls[2])  # 1-based to 0-based
    tmp=[]
    for i in range(args.refseq_start, args.refseq_end):
        tmp.append(d[i])
    mean_depth=np.mean(tmp)
    for i in range(args.refseq_start, args.refseq_end - params.k + 1):
        if i % params.slide_bin == 0:
            tmp=[]
            for apos in range(params.k):
                tmp.append(d[i + apos])
            tmp_depth=np.mean(tmp)
            coeff[i]= mean_depth / tmp_depth


def sam_to_kmer(args, params, filenames, coeff):
    log.logger.debug('started')
    try:
        def complement(string):
            seq_c=string.translate(str.maketrans('ATGC', 'TACG'))[::-1]
            return seq_c
        
        def count_clip(cigar):
            length=0
            tmp=''
            for c in cigar:
                if not c in cigar_op:
                    tmp += c
                elif c == 'S':
                    length += int(tmp)
                    tmp=''
                else:
                    tmp=''
            return length
        
        def calc_ref_len(cigar):
            length=0
            tmp=''
            for c in cigar:
                if not c in cigar_op:
                    tmp += c
                elif c in cigar_ref_retain:
                    length += int(tmp)
                    tmp=''
                else:
                    tmp=''
            return length
        
        def parse_seq_to_kmers(seq, cigar, genome_start, genome_end, chr):
            # parse cigar
            parsed_cigar=[]
            tmp=''
            for c in cigar:
                if not c in cigar_op:
                    tmp += c
                else:
                    parsed_cigar.append([c, int(tmp)])
                    tmp=''
            # clip seq
            left,right=0,0
            if 'S' in cigar:
                if parsed_cigar[0][0] == 'S':
                    left=parsed_cigar[0][1]
                if parsed_cigar[-1][0] == 'S':
                    right=parsed_cigar[-1][1]
            ref_pos=genome_start
            read_pos=left
            pos_d={}
            for c,l in parsed_cigar:
                if c == 'M':
                    for _ in range(l):
                        pos_d[ref_pos]=read_pos
                        read_pos += 1
                        ref_pos += 1
                elif c == 'I':
                    for _ in range(l):
                        read_pos += 1
                elif c == 'D':
                    for _ in range(l):
                        pos_d[ref_pos]=read_pos
                        ref_pos += 1
            pos_d[ref_pos]= read_pos
            tmp=[]
            if (len(seq) - (left + right)) >= params.k:
                for i in range(genome_start, genome_end - params.k + 1):
                    if i % params.slide_bin == 0:
                        kmer_seq=seq[pos_d[i]:pos_d[i + params.k]]
                        ref='%s:%d-%d' % (chr, i, i + params.k)  # 0-based start; 1-based end
                        tmp.append([kmer_seq, ref, coeff[i]])
            return tmp

        import pysam
        if args.c is None:
            infile=pysam.AlignmentFile(args.b, 'rb')
        else:
            infile=pysam.AlignmentFile(args.c, 'rc', reference_filename=args.fa)
        
        kmers_set=set()
        kmers_ls=[]
        for line in infile:
            line=line.tostring()
            ls=line.strip().split('\t')
            if ls[2] == args.refseq_id:   # only specified refseq
                ref_len=calc_ref_len(ls[5])
                map_start= int(ls[3]) - 1  # 0-based
                map_end= map_start + ref_len
                if args.refseq_start <= map_start and map_end <= args.refseq_end:
                    if int(ls[1]) < 256:  # only use primary alignment, discard supplementary and duplicate
                        b=bin(int(ls[1]))
                        if b[-2] == '1':   # properly paired
                            for i in ls[::-1]:
                                if 'NM:i:' in i:
                                    nm= int(i.replace('NM:i:', ''))
                            if 'S' in ls[5]:
                                soft_clip_len=count_clip(ls[5])
                            else:
                                soft_clip_len=0
                            if nm <= params.max_mut and soft_clip_len <= params.max_clip_len:
                                kmers_l=parse_seq_to_kmers(ls[9], ls[5], map_start, map_end, ls[2])  # 0-based
                                kmers_ls.extend(kmers_l)
                                for kmer,_,_ in kmers_l:
                                    kmers_set.add(kmer)
        kmers_set=sorted(list(kmers_set))
        kmers_d={}
        counts_d={}
        for kmer in kmers_set:
            kmers_d[kmer]=[]
            counts_d[kmer]=[]
        for kmer,pos,norm_count in kmers_ls:
            kmers_d[kmer].append(pos)
            counts_d[kmer].append(norm_count)
        out=[]
        for kmer in kmers_d:
            pos_set=sorted(list(set(kmers_d[kmer])))
            norm_counts=sum(counts_d[kmer])
            attr=[]
            for pos in pos_set:
                attr.append('%s=%d' % (pos, kmers_d[kmer].count(pos)))
            out.append('%s\t%d\t%f\t%s\n' % (kmer, len(kmers_d[kmer]), norm_counts, ';'.join(attr)))
        with open(filenames.summary, 'w') as outfile:
            outfile.write(''.join(out))
            
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
