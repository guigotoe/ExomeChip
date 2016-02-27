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
    inf = sys.argv[1]
    out = open('duplicates.txt','w')

    for i in open(inf,'r'):
        if re.findall('Duplicate',i): out.write(i.split("[")[1].strip(" ]\n")+"\n")

if __name__ == '__main__':main()