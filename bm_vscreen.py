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

#from Bio import pairwise2
from ktio.ktio import readFastq

VectorContamination = namedtuple('VectorContamination', 'start end location coverage'.split(' '))
BlastHSP = namedtuple('BlastHSP', 'query subject pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen positive gaps ppos frames staxids salltitles sstrand qseq sseq'.split(' '))
FQRead = namedtuple('FQRead', 'id seq qual'.split(' '))
BACAssignedRead = namedtuple('BACAssignedRead', 'id seq qual bacid'.split(' '))
FQSplitRead = namedtuple('FQSplitRead', 'id seqs quals bacid'.split(' '))
NamedFile = namedtuple('NamedFile', 'filename file'.split(' '))

BLASTCMD = '{} -perc_identity 95 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -db {} -outfmt "6 std qlen slen positive gaps ppos frames staxids salltitles sstrand qseq sseq"'
ECODB = './eco.fa'
CONTAMINANTS = {'NC_010473.1'}
VECTOR = 'pIndigoBeloHindIII'

def mask(seq, start, end, soft=False):
    return ''.join((seq[:start - 1], seq[start - 1:end].lower() if soft else (end - start + 1) * 'N', seq[end:]))

def runBlastAll(query, blastcmd, blastdb, blaster):
    pr = sub.Popen(blastcmd.format(blaster, blastdb), shell=True, stdin=sub.PIPE, stderr=sub.PIPE, stdout=sub.PIPE)
    out, err = pr.communicate(query.encode())
    out = out.decode().strip()
    for hsp in out.split('\n'):
       if hsp:
           yield BlastHSP(*(hsp.split('\t')))

def vscreen(bacid, fn):
    print('VSCREEN:', bacid, fn, file=sys.stderr)
    reads = dict((_id[1:], [_seq, _qual, _seq]) for _id,_seq,_qual in readFastq(fn))

    '''
    all_vc, vc, query = dict(), dict(), str()
    novc = set()
    while True:
        print('VC', bacid, list(vc.keys())[:10], file=sys.stderr)
        for k in vc:
            all_vc.setdefault(k, list()).append(vc[k])
            _seq = mask(reads[k][0], vc[k])
            reads[k] = (_seq, reads[k][1])
        contaminated = set(reads).intersection(vc)
        if not reads or (query and not contaminated):
            break
        novc = set(reads).difference(all_vc) if all_vc else set()
    '''


    vloc, contaminated = dict(), set()
    dbHits, query = set(), str()
    vc, nohits = set(), set()
    allhits = set()
    first = True
    while True:
        allhits.update(dbHits)
        if not first and not dbHits:
            break
        first = False
        nohits = set(reads).difference(allhits) if allhits else set()
        dbHits = set()
        query = ''.join(('>{}\n{}\n'.format(_id, reads[_id][2]) for _id in reads if _id not in nohits))
        for hsp in runBlastAll(query, BLASTCMD, ECODB, 'blastn'):
            alen = int(hsp.qend) - int(hsp.qstart) + 1
            if alen < 200:
                continue
            dbHits.add(hsp.query)
            reads[hsp.query][2] = mask(reads[hsp.query][2], int(hsp.qstart), int(hsp.qend))
            if hsp.query in CONTAMINANTS and alen/qlen >= 0.9:
                contaminated.add(hsp.query)
            else:
                vloc.setdefault(hsp.query, list()).append((int(hsp.qstart), int(hsp.qend)))

    with open(fn.replace('.fq', '.masked.fq'), 'wt') as masked_out, open(fn.replace('.fq', '.clean.fq'), 'wt') as clean_out:
        for _id in sorted(reads):
            cov = reads[_id][2].count('N') / len(reads[_id][2])
            if _id in contaminated or cov >= 0.9:
                print(_id, "CONTAMINATED", cov, file=sys.stderr, sep='\t')

            else:
                seq, _, masked = map(list, reads[_id])
                for start, end in vloc.get(_id, list()):
                    masked[start-1:end] = list(map(str.lower, seq[start-1:end]))
                reads[_id][2] = ''.join(masked)
                print('@{}\n{}\n+\n{}'.format(_id, reads[_id][2], reads[_id][1]), file=masked_out)
                unmasked = [match for match in sorted(list(re.finditer('[acgtACGT]+', reads[_id][2])), key=lambda x:len(x.group()), reverse=True) if len(match.group())>=200][:2]
                if unmasked:
                    outid = _id.replace('/css', '.1/css')
                    print('@{}\n{}\n+\n{}'.format(outid, unmasked[0].group(), reads[_id][1][unmasked[0].start():unmasked[0].end()]), file=clean_out)
                    if len(unmasked) == 2:
                        outid = _id.replace('/css', '.2/css')
                        print('@{}\n{}\n+\n{}'.format(outid, unmasked[1].group(), reads[_id][1][unmasked[1].start():unmasked[1].end()]), file=clean_out)
    return [bacid + ' done']

def processFileMP3(outfiles, nthreads):
    pool = mp.Pool(processes=nthreads)
    results = [pool.apply_async(vscreen, args=(bacid, outfiles[bacid])) for bacid in outfiles]
    for query_results in (p.get() for p in results):
        pass

if __name__ == '__main__':

    INPUTFILES = dict((os.path.basename(f).strip('.reads.fq'), f)
                      for f in glob.glob(os.path.join(sys.argv[1], '*.reads.fq')))
    processFileMP3(INPUTFILES, 3)
