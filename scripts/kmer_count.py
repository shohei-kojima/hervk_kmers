#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys
import log,traceback


cigar_op={'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
cigar_ref_retain={'M', 'D', 'N', '=', 'X'}
cigar_read_retain={'M', 'I', '=', 'X'}


def sam_to_kmer(args, params, filenames):
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
        
        def parse_seq_to_kmers(seq, cigar, genome_start, chr):
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
                        pos_d[read_pos]=ref_pos
                        read_pos += 1
                        ref_pos += 1
                elif c == 'I':
                    for _ in range(l):
                        pos_d[read_pos]=ref_pos
                        read_pos += 1
                elif c == 'D':
                    for _ in range(l):
                        ref_pos += 1
            pos_d[read_pos]= ref_pos
            tmp=[]
            if (len(seq) - (left + right)) >= params.k:
                for i in range(left, len(seq) - right - params.k + 1, params.slide_bin):
                    kmer_seq=seq[i:i + params.k]
                    ref='%s:%d-%d' % (chr, pos_d[i], pos_d[i + params.k])  # 0-based start; 1-based end
                    tmp.append([kmer_seq, ref])
            return tmp

        import pysam
        infile=pysam.AlignmentFile(args.b, 'rb')
        
        kmers_set=set()
        kmers_ls=[]
        for line in infile:
            line=line.tostring()
            ls=line.strip().split('\t')
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
                        if b[-5] == '1':
                            seq=complement(ls[9])
                        else:
                            seq=ls[9]
                        kmers_l=parse_seq_to_kmers(seq, ls[5], int(ls[3]) - 1, ls[2])  # 0-based
                        kmers_ls.extend(kmers_l)
                        for kmer,_ in kmers_l:
                            kmers_set.add(kmer)
        kmers_set=sorted(list(kmers_set))
        kmers_d={}
        for kmer in kmers_set:
            kmers_d[kmer]=[]
        for kmer,pos in kmers_ls:
            kmers_d[kmer].append(pos)
        out=[]
        for kmer in kmers_d:
            pos_set=sorted(list(set(kmers_d[kmer])))
            attr=[]
            for pos in pos_set:
                attr.append('%s=%d' % (pos, kmers_d[kmer].count(pos)))
            out.append('%s\t%d\t%s\n' % (kmer, len(kmers_d[kmer]), ';'.join(attr)))
        with open(filenames.summary, 'w') as outfile:
            outfile.write(''.join(out))
            
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
