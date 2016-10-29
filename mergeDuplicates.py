#!/usr/bin/env python

import os, sys, distutils
import argparse

parser=argparse.ArgumentParser(description='This script takes in a tab-delimited file or STDIN (use "-"), '
                                           'and a column number and merges lines that show the same entry on '
                                           'that column number. It searches and merges if the other entries are the '
                                           'same, otherwise, it uses a ";" to merge',
                               usage='mergeDuplicates.py <file> <colnum>',
                               epilog="EXAMPLE: mergeDuplicates.py test.txt 3 > testout.txt")
parser.add_argument('file', nargs='+', help='tab-delimited input file')
parser.add_argument('colnum', nargs=1, default=1, help='integer; 1-based')
parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='STDOUT; Note that the '
                                                                                               'order '
                                                                                               'is NOT preserved; '
                                                                                               'sorted in '
                                                                                               'alphanumeric order by '
                                                                                               'desired column')
## help if no arguments or -h/--help
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args=parser.parse_args()

if __name__ == '__main__':

    ## variables
    if sys.argv[1] == '-':
        lines = sys.stdin.readlines()
    else:
        f1 = open(sys.argv[1])
        lines = f1.xreadlines()

    colnum = int(sys.argv[2]) - 1
    myuniqlines = {}

    for line in lines:
        fields = line.rstrip().split('\t')

        ## if there is a duplicate row (based on the column),
        ## merge such that anything that's not duplicated, separate by ';'
        ## if everything is duplicated, it's essentially a sort and uniq
        if fields[colnum] in myuniqlines:
            for i in range(len(fields)):
                if myuniqlines[ fields[colnum] ].rstrip().split('\t')[i] != fields[i]:
                    mfields = myuniqlines[ fields[colnum] ].rstrip().split('\t')
                    mfields[i] = mfields[i] + ';' + fields[i]
                    newline = '\t'.join(mfields)
                    myuniqlines[fields[colnum]] = newline

        else:
            myuniqlines[ fields[colnum] ] = line.rstrip()


    ## print stored lines
    for key in sorted(myuniqlines):
        print myuniqlines[key]




