#!/usr/bin/python

######################################################
## By Guillermo Torres PhD. - IKMB                 ##
####################################################
#
#* last update 5.2015
#
#** Pendent for next update:
#**     - check if there are repeated markers
#**     - Graph of New markers in function of the chromosome where they are located
#**
#******************

## Libraries importing

import pandas as pd
import sys, os
from optparse import OptionParser
import matplotlib.pyplot as plt

## Parsing Options ##

parser = OptionParser()
usage = """
This scripts is able to check two MAP files from genotype data (GWAS chip experiment) and exclude a SNPs from further analysis.
It output is a marker differences log file as well as the not shared markers exclusion (it has be specified
by the options).
* Please ensure that you have a backup map file before change to a new one!
REQUIREMENTS:\n** pandas and matplotlib python libraries
\t%prog --if1 PATH/first.map --if2 PATH/second.map [args] \n\t%prog {-h --help}\n
This is part of ExomeChip analysis by Guillermo G. Torres (g.torres@ikmb.uni-kiel.de) - Genetics of Human Longevity
Group - IKMB - Christian-Albrechts-Universitat zu Kiel
"""
parser = OptionParser(usage=usage,version="%prog 0.1")
parser.add_option("-1","--if1",type="string",action="store",dest="map1",#default='./1.map',
                      help="Map file One - Mandatory")
parser.add_option("-2","--if2",type="string",action="store",dest="map2",#default='./2.map',
                      help="Map file Two - Mandatory")
parser.add_option("-o","--out",type="string",action="store",dest="outlog",default='dataExpl.log',
                      help="Out logfile name - default name = dataExpl.log")
parser.add_option("-c", action="store_true",dest="c",default=True,
                      help = "Change original map file content (exclude SNPs) - default True")
parser.add_option("-l", type="string",action="store",dest="lb",default='cases,controls',
                      help = "File labels, put them separated with comma - ie. a,b")
parser.add_option("-v", action="store_true",dest="verbose",default=False,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()

if len(args)!=0 or o.map1 == None or o.map2 == None:    # assertion if args. were put well.
    parser = OptionParser()
    parser.error("""Incorrect number of arguments\n\
\t%prog --if1 PATH/first.map --if2 PATH/second.map [args] \n\t%prog {-h --help}\n
This is part of ExomeChip analysis by Guillermo G. Torres (g.torres@ikmb.uni-kiel.de) - Genetics of Human Longevity
Group - IKMB - Christian-Albrechts-Universitat zu Kiel""")
else:   # assertion if files exist
    if os.path.isfile(o.map1) is False and os.path.isfile(o.map2) is False: parser.exit("Error:\n\tFiles 1 and 2 do not exist: %s %s\n"%(o.map1,o.map2))
    elif os.path.isfile(o.map1) is False: parser.exit("Error:\n\tFile 1 does not exist: %s\n"%o.map1)
    elif os.path.isfile(o.map2) is False: parser.exit("Error:\n\tFile 2 does not exist: %s\n"%o.map2)
if o.lb != None:
    if len(o.lb.split(',')) == 2: o.lb = o.lb.split(',');pass
    elif len(o.lb.split(',')) == 1: sys.stdout.write("Warning!:\n\tJust 1 label was provided: %s\n"%(o.lb.split(',')[0])); o.lb = [o.lb.split(',')[0],'map_2']
    elif len(o.lb.split(',')) > 2: sys.stdout.write("Warning!:\n\tMore than 2 labels were provided: %s\n"%(' '.join(o.lb.split(','))));o.lb = o.lb.split(',')
else:
    o.lb='map_1,map_2'
## Global variables ##
log_out = open(o.outlog, 'w')
log_out.close()     # Re-starting the log file...
if os.path.dirname(o.map1) == '.' or os.path.dirname(o.map1) == '':pm1=os.getcwd()+'/'
else:pm1=os.path.dirname(o.map1)+'/'
if os.path.dirname(o.map2) == '.' or os.path.dirname(o.map2) == '':pm2=os.getcwd()+'/'
else:pm2=os.path.dirname(o.map2)+'/'

## Main ##

def main():
    log_out = open(o.outlog, 'a')
    m1 = pd.read_table(o.map1,names=['Chr','SNP_id','GDis','Pos']).fillna(0)
    m2 = pd.read_table(o.map2,names=['Chr','SNP_id','GDis','Pos']).fillna(0)
    df = pd.concat([m1,m2],keys=o.lb[:2])
    df['Chr'] = df['Chr'].astype(str);df['Pos'] = df['Pos'].astype(int)
# First we grouping the similar ones... each group will be append two
    df_gpby = df.groupby(list(df.columns)) # same as : df_gpby = df.groupby(by=['Chr','SNP','GDis','Pos'])

#*** sharing ***#
    idy = []
    for y in df_gpby.groups.values():
        if len(y) == 2:idy.extend([y[0],y[1]])
    shared = df.reindex(idy).sort_index()

#*** all diferences ***#
    idx = [x[0] for x in df_gpby.groups.values() if len(x) == 1]    # appending the indexes from the different ones in idx array...
    diffs = df.reindex(idx).sort_index()
    log_out.write("The following summarizes all differences between the two map files:\n\n")    # All
    diffs.to_csv(log_out,sep="\t",index_label=['mapfile','idx'])

#*** diferences by: ***#
    each = {'Pos':['Chr','SNP_id','GDis'],'SNP_id':['Chr','Pos','GDis'],'Chr':['SNP_id','Pos','GDis'],'GDis':['SNP_id','Pos','Chr']}
    each2 = {'Chr & SNP_id':'Pos','Chr & Pos':'SNP_id'}
    news = diffs
    stats = {'Shared':len(shared)/2}
    for v in each.values():
        gpby = news.groupby(by=v)
        key =each.keys()[each.values().index(v)]
        gidx = []
        for i in gpby.groups.values():
            if len(i) == 2:
                if i[0][0] != i[1][0]:gidx.extend([i[0],i[1]])
        dby = news.reindex(gidx)
        if len(dby) > 0:
            log_out.write("\nDifferences only by: %s\n"%key)
            dby.to_csv(log_out,sep="\t",index_label=['mapfile','idx'])
            stats[key]=len(dby)/2
        else:stats[key]=0
        news=news.drop(gidx)

    for v in each2.values():
        gpby = news.groupby(by=v)
        key =each2.keys()[each2.values().index(v)]
        gidx = []
        for i in gpby.groups.values():
            if len(i) == 2:gidx.extend([i[0],i[1]])
        dby = news.reindex(gidx)
        if len(dby) > 0:
            log_out.write("\nDifferences both by: %s\n"%key)
            dby.to_csv(log_out,sep="\t",index_label=['mapfile','idx'])
            stats[key]=len(dby)/2
        else:stats[key]=0
        news=news.drop(gidx)
    if len(news) > 0:
        log_out.write("\nNew Markers!:\n")
        news.to_csv(log_out,sep="\t",index_label=['mapfile','idx'])
        stats['New']=len(news)
        stats["New-%s"%o.lb[0]]=len(news.xs(o.lb[0]))
        stats["New-%s"%o.lb[1]]=len(news.xs(o.lb[1]))
        c={}
        for i in news.groupby('Chr').groups.values():c[news.ix[i,'Chr'][0]]=len(i)
        chr = pd.DataFrame.from_dict(c, orient="index")
        chr = chr.rename(columns = {0:'New-Marker'}).sort_index()
        chr.plot(kind='bar')
        plt.savefig('NewMarkers.jpeg',dpi=150,size=(6,6),bbox_inches='tight')
    else: stats['New'] = 0

#*** Ploting results!! ***#
    sdf = pd.DataFrame.from_dict(stats, orient="index")
    sdf=sdf.rename(columns = {0:'Count'}).sort_index()
    log_out.write("\n");sdf.to_csv(log_out,sep="\t",index_label=['Diferences'])
    sdf.plot(kind='bar')
    plt.savefig('stats.jpeg',dpi=150,size=(6,6),bbox_inches='tight')

    #plt.show()
    log_out.close()

#** Writing **#
    if o.c is True:
        for i in idx:df.ix[i,'Pos']=-df.ix[i,'Pos']     # Changing the Position to negative numbers
        os.rename(pm1+os.path.basename(o.map1), pm1+'deBckup_'+os.path.basename(o.map1))
        nfm1 = open(pm1+os.path.basename(o.map1), 'w')
        df.ix[o.lb[0]].to_csv(nfm1,sep="\t",index=False,header=False)
        os.rename(pm2+os.path.basename(o.map2), pm2+'deBckup_'+os.path.basename(o.map2))
        nfm2 = open(pm2+os.path.basename(o.map2), 'w')
        df.ix[o.lb[1]].to_csv(nfm2,sep="\t",index=False,header=False)
        print os.path.basename(o.map1), os.path.basename(o.map2)
        nfm1.close();nfm2.close()

    print('Done!')

## Functions ##
def report_diff(x):
    print x
    #return x[0] if x[0] == x[1] else '{} | {}'.format(*x)

if __name__ == '__main__':main()