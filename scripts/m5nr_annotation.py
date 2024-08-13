#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get the annotations from a list of md5s.
Miguel Romero github.com/romeromig
"""

import pickle
import sys
import os

db_dict = {
        'm5nr_genbank.dict':['GenBank','97'],
        'm5nr_refseq.dict':['RefSeq','60'],
        'm5nr_seed.dict':['SEED','20'],
        'm5nr_kegg.dict':['KEGG','6']
        }

if len(sys.argv[1:]) == 0:
    print("""
                Get the annotations from a list of md5s.
                usage: python3 m5nr_annotation.py database.dict md5s.txt output
                Available databases: 
                m5nr_genbank.dict = GenBank
                m5nr_refseq.dict = RefSeq
                m5nr_seed.dict = SEED
                m5nr_kegg.dict = KEGG GENES
                Miguel F. Romero github.com/miferg
                """
                )
    sys.exit()

argv = sys.argv[1:]



print()
print('''
Get the annotations from a list of md5s.
Miguel F. Romero github.com/miferg
''')

dbfile = os.path.abspath(argv[0]).split('/')[-1]

if dbfile in db_dict.keys():
    print('Database: '+ db_dict[dbfile][0])
else:
    print('Please use one of the available databases.\n')
    sys.exit()

print('Loading database, this may take a few minutes (and ~'+ db_dict[dbfile][1] +' GBs of memory).')

database = pickle.load(open(os.path.abspath(argv[0]), "rb"))

print('Loading input.')

handle = open(argv[1],'r')
md5s = handle.read().split('\n')[:-1]
handle.close()

print('Writing output in '+ argv[2])
outfile = open(argv[2], 'w')

absent = []

if dbfile == 'm5nr_genbank.dict':
    for md5 in md5s:
        try:
            outfile.write(md5 +'\t'+ database[md5] +'\n')
        except:
            absent.append(md5)
    outfile.close()

elif dbfile == 'm5nr_refseq.dict':
    for md5 in md5s:
        try:
            outfile.write(md5 +'\t'+ '\t'.join(list(database[md5])) +'\n')
        except:
            absent.append(md5)
    outfile.close()

elif dbfile == 'm5nr_seed.dict':
    for md5 in md5s:
        try:
            ssids = database[md5]
            for line in ssids:
                outfile.write(md5 +'\t'+ line)
        except:
            absent.append(md5)
    outfile.close()

elif dbfile == 'm5nr_kegg.dict':
    for md5 in md5s:
        try:
            outfile.write(md5 +'\t'+ database[md5] +'\n')
        except:
            absent.append(md5)
    outfile.close()

if len(absent) > 0:
    print('There were '+ str(len(absent)) +' md5s missing in the annotation database.')
    print('Check the m5nr_annotation_missing.txt file.')
    outfile = open('m5nr_annotation_missing.txt','w')
    outfile.write('\n'.join(absent) +'\n')
    outfile.close()

print('Done!\n')
