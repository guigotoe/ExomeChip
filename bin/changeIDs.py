#!/usr/bin/python

######################################################
## By Guillermo Torres PhD. - IKMB                 ##
####################################################
#
#* last update 6.2015
#******************

## Libraries importing

import sys, re

def main():
    inf = "/home/torres/Documents/Projects/ExomeChip/bin/ExonChipProcessing/extra_files/chr23_26.txt"
    inf2 = "/home/torres/Documents/Projects/ExomeChip/bin/ExonChipProcessing/extra_files/snpid_updated.txt"
    outf = "/home/torres/Documents/Projects/ExomeChip/bin/ExonChipProcessing/extra_files/chr23_26_newIDs.txt"
    out = open(outf,'w')
    old = {}
    new = {}
    for i in open(inf2,'r'):
        i = i.strip("\n").split("\t")
        old[i[0]] = i[1]
    for j in open(inf,'r'):
        j = j.strip("\n\r")
        if j in old: out.write("%s\n"%old[j])
        else:pass
    print ("%s successfully written !"%outf)
if __name__ == '__main__':main()