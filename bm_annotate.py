#!/usr/bin/env python
# coding=utf-8
import os
import sys
import re
import csv
import glob
import subprocess as sub
import multiprocessing as mp

from collections import Counter, namedtuple

from ktio.ktio import readFastq

BLASTCMD = ' {} -query {} -perc_identity 98 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -db {} -outfmt "6 std qlen slen positive gaps ppos frames staxids salltitles sstrand qseq sseq"'
BlastHSP = namedtuple('BlastHSP', 'query subject pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen positive gaps ppos frames staxids salltitles sstrand qseq sseq'.split(' '))

VECTOR = 'pIndigoBeloHindIII'

def runBlastAll(query, blastcmd, blastdb, blaster):
    pr = sub.Popen(blastcmd.format(blaster, query, blastdb), shell=True, stdin=sub.PIPE, stderr=sub.PIPE, stdout=sub.PIPE)
    out, err = pr.communicate()
    out = out.decode().strip()
    # print(out, err)
    for hsp in out.split('\n'):
       if hsp:
           yield BlastHSP(*(hsp.split('\t')))



def annotate(fn, db, feature):
    # seqs = dict((_id[1:], [_seq, _qual, _seq]) for _id,_seq,_qual in readFastq(fn))
    for hsp in runBlastAll(fn, BLASTCMD, db, 'blastn'):
        #print(hsp)
        qstart, qend = int(hsp.qstart), int(hsp.qend)
        sstart, send = sorted(map(int, (hsp.sstart, hsp.send)))
        alen = qend - qstart + 1
        if alen >= 200:

            attr = 'id={};slen={};qlen={};scov={:.3f}'.format(hsp.subject, hsp.slen, hsp.qlen, (send - sstart + 1)/int(hsp.slen))
            print(hsp.query, 'qad-ann', feature, qstart, qend, hsp.pident, '+' if hsp.sstrand == 'plus' else '-', '.', attr, sep='\t')


if __name__ == '__main__':
    annotate(sys.argv[1], sys.argv[2], sys.argv[3])
