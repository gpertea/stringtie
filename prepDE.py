#!/usr/bin/env python
import re, csv, sys, os, glob, warnings, itertools
from math import ceil
from optparse import OptionParser
from operator import itemgetter
#note that the gtf files in the sample folders have same # of lines, just different order(?)

parser=OptionParser(description='Generates two CSV files containing the count matrices for genes and transcripts, using the coverage values found in the output of `stringtie -e`')
parser.add_option('-i', '--input', '--in', default='ballgown', help="the parent directory of the sample sub-directories [default: %default]")
parser.add_option('-g', default='gene_count_matrix.csv', help="where to output the gene count matrix [default: %default")
parser.add_option('-t', default='transcript_count_matrix.csv', help="where to output the transcript count matrix [default: %default]")
parser.add_option('-l', '--length', default=75, type='int', help="the average read length [default: %default]")
parser.add_option('-p', '--pattern', default=".", help="a regular expression that selects the sample subdirectories")
parser.add_option('-c', '--cluster', action="store_true", help="whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)")
parser.add_option('-s', '--string', default="MSTRG", help="if a different prefix is used for geneIDs assigned by StringTie [default: %default]")
parser.add_option('-k', '--key', default="prepG", help="if clustering, what prefix to use for geneIDs assigned by this script [default: %default]")
parser.add_option('--legend', default="legend.csv", help="if clustering, where to output the legend file mapping transcripts to assigned geneIDs [defaukt: %default]")
(opts, args)=parser.parse_args()
#opts.cluster=True
if not os.path.isdir(opts.input):
  parser.print_help()
  print " "
  print "Error: sub-directory '%s' not found!" % ( opts.input)
  sys.exit(1)
RE_GENE_ID=re.compile('gene_id "([^"]+)"')
RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')
RE_COVERAGE=re.compile('cov "([\-\+\d\.]+)"')
RE_STRING=re.compile(re.escape(opts.string))

samples = [i for i in next(os.walk(opts.input))[1] if re.search(opts.pattern,i)]
if len(samples) == 0:
  parser.print_help()
  print " "
  print "Error: no GTF files found under ./%s !" % ( opts.input)
  sys.exit(1)

samples.sort()


def is_transcript(x):
  return len(x)>2 and x[2]=="transcript"

def getGeneID(s, ctg, tid):
  r=RE_GENE_ID.search(s)
  if r: return r.group(1)
  r=RE_GENE_NAME.search(s)
  if r: return ctg+'|'+r.group(1)
  return tid

def getCov(s):
  r=RE_COVERAGE.search(s)
  if r: 
    v=float(r.group(1))
    if v<0.0: v=0.0
    return v
  return 0.0

def is_overlap(x,y): #NEEDS TO BE INTS!
  return x[0]<=y[1] and y[0]<=x[1]


def t_overlap(t1, t2): #from badGenes: chromosome, strand, cluster, start, end, (e1start, e1end)...
    if t1[0] != t2[0] or t1[1] != t2[1] or t1[5]<t2[4]: return False
    for i in xrange(6, len(t1)):
        for j in xrange(6, len(t2)):
            if is_overlap(t1[i], t2[j]): return True
    return False


read_len=opts.length
t_count_matrix, g_count_matrix=[],[]

##Get ready for clustering, stuff is once for all samples##
geneIDs={} #key=transcript, value=cluster/gene_id

for s in samples:
    badGenes=[] #list of bad genes (just ones that aren't MSTRG)
    
    try:
        with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f:
            split=[l.split('\t') for l in f.readlines()]
        for i,v in enumerate(split):
            if is_transcript(v):
                t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
                try:
                  g_id=getGeneID(v[len(v)-1], v[0], t_id)
                except:
                  print "Problem at line:\n:%s\n" % (v)
                  print "i='%s', len(v)=%s" % (i, len(v));
                  sys.exit(1)
                geneIDs.setdefault(t_id, g_id)
                if not RE_STRING.match(g_id):
                    badGenes.append([v[0],v[6], t_id, g_id, min(int(v[3]),int(v[4])), max(int(v[3]),int(v[4]))]) #chromosome, strand, cluster/transcript id, start, end
                    j=i+1
                    while j<len(split) and split[j][2]=="exon":
                        badGenes[len(badGenes)-1].append((min(int(split[j][3]), int(split[j][4])), max(int(split[j][3]), int(split[j][4]))))
                        j+=1

    except StopIteration:
        warnings.warn("Didn't get a GTF in that directory. Looking in another...")

    else:
        break

##THE CLUSTERING BEGINS!##
if opts.cluster and len(badGenes)>0:
    clusters=[] #lists of lists (could be sets) or something of transcripts
    badGenes.sort(key=itemgetter(3)) #sort by start coord...?
    i=0
    while i<len(badGenes): #rather un-pythonic
        temp_cluster=[badGenes[i]]

        k=0
        while k<len(temp_cluster):
            j=i+1
            while j<len(badGenes):
                if t_overlap(temp_cluster[k], badGenes[j]):
                    temp_cluster.append(badGenes[j])
                    del badGenes[j]
                else:
                    j+=1
            k+=1            
        if len(temp_cluster)>1:
            clusters.append([t[2] for t in temp_cluster])
        i+=1

    print len(clusters)

    sortedClusters=clusters
    for c in sortedClusters:
        c.sort()

    sortedClusters.sort(key=itemgetter(0))
    print sortedClusters
    legend=[]
    for u,c in enumerate(clusters):
        my_ID=opts.key+str((u+1))
        print my_ID
        legend.append(list(itertools.chain.from_iterable([[my_ID],c]))) #my_ID, clustered transcript IDs
#        print c
        for t in c:
            geneIDs[t]=my_ID            
##            geneIDs[t]="|".join(c) #duct-tape transcript IDs together, disregarding ref_gene_names and things like that

    print legend
    
    with open(opts.legend, 'w') as l_file:
        my_writer=csv.writer(l_file)
        my_writer.writerows(sortedClusters)
    
geneDict={} #key=gene/cluster, value=dictionary with key=sample, value=summed counts
t_dict={}
for q, s in enumerate(samples):
    print q,s

    try:
        with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f: #grabs first .gtf file it finds inside the sample subdirectory

##        split=[t[:len(t)-1]+t[len(t)-1].split(";") for t in split]
##        split=[t[:len(t)-1] for t in split] #eliminate '\n' at end
##        split=[[e.lstrip() for e in t] for t in split]
        #should consider making stuff into dictionaries, maybe each split line

##            transcriptList=[]
            transcript_len=0
            for l in f:
                if l.startswith("#"):
                    continue
                v=l.split('\t')
                if v[2]=="transcript":
                    if transcript_len>0:
##                        transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
                        t_dict.setdefault(t_id, {})
                        t_dict[t_id].setdefault(s, int(ceil(coverage*transcript_len/read_len)))
                    t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
                    #g_id=RE_GENE_ID.search(v[len(v)-1]).group(1)
                    g_id=getGeneID(v[len(v)-1], v[0], t_id)
                    #coverage=float(RE_COVERAGE.search(v[len(v)-1]).group(1))
                    coverage=getCov(v[len(v)-1])
                    transcript_len=0
                if v[2]=="exon":
                    transcript_len+=int(v[4])-int(v[3])

##            transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
            t_dict.setdefault(t_id, {})
            t_dict[t_id].setdefault(s, int(ceil(coverage*transcript_len/read_len)))

    except StopIteration:
        warnings.warn("No GTF file found in "+os.path.join(opts.input,s))

##        transcriptList.sort(key=lambda bla: bla[1]) #gene_id

    for i,v in t_dict.iteritems():
##        print i,v
        geneDict.setdefault(geneIDs[i],{}) #gene_id
        geneDict[geneIDs[i]].setdefault(s,0)
        geneDict[geneIDs[i]][s]+=v[s]

with open(opts.t, 'w') as csvfile:
    my_writer=csv.DictWriter(csvfile, fieldnames=[""]+samples)
    my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
    for i in t_dict:
        t_dict[i][""]=i
        my_writer.writerow(t_dict[i])

with open(opts.g, 'w') as csvfile:
    my_writer=csv.DictWriter(csvfile, fieldnames=[""]+samples)
##    my_writer.writerow([""]+samples)
##    my_writer.writerows(geneDict)
    my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
    for i in geneDict:
        geneDict[i][""]=i #add gene_id to row
        my_writer.writerow(geneDict[i])
