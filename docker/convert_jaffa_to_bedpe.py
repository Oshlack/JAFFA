#!/usr/bin/env python
# Author: Rebecca Evans (rebecca.louise.evans@gmail.com)
import sys
import string

with open(sys.argv[1], 'r') as jaffaOutput:
    next(jaffaOutput) # skip header
    for line in jaffaOutput:
        vals = line.strip().split(',')
        chrL = vals[2].strip('"')
        posL = vals[3].strip('"')
        strandL = vals[4].strip('"')
        geneL = vals[1].split(':')[0].strip('"')
        chrR = vals[5].strip('"')
        posR = vals[6].strip('"')
        strandR = vals[7].strip('"')
        geneR = vals[1].split(':')[1].strip('"')
        bedpe = '\t'.join([chrL, posL, str(int(posL)+1), chrR, posR,
            str(int(posR)+1), '>>'.join([geneL, geneR]), '0', strandL, strandR])
        print(bedpe)
