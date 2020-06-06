#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging

'''
time python main.py -overwrite -b /home/kooojiii/results/2020/smrv/200531_1/tmp.bam
'''


# version
version='2020/06/02'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')
parser.add_argument('-fa', metavar='str', type=str, help='When you use -c option, specify reference genome which are used when input reads were mapped. Example: hg38.fa')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-v', '--version', help='Print version.', action='store_true')
args=parser.parse_args()
args.version=version


# start
import init
init.init(args, version)


# logging
import log
args.logfilename='for_debug.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
log.logger.debug('This is version %s' % version)
print()
log.logger.info('Initial check started.')
initial_check.check(args, sys.argv)


# set up
import setup
setup.setup(args, init.base)
params=setup.params


# output file names
import utils
filenames=utils.empclass()
filenames.summary=os.path.join(args.outdir, 'kmer_counts.txt')


# 1. K-mer count
import kmer_count
log.logger.info('K-mer counting started.')
kmer_count.sam_to_kmer(args, params, filenames)

log.logger.info('All analysis finished!')
