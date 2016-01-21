#!/usr/bin/env python

# plot_tripseq_transcript.py 

# plot the indicated transcript or gene based input polysome sequencing (TrIPseq) data

# Intended to be applied to fractionated polysome sequencing data (a la tripseq) but theoretically can work for any data involving transcript-specific distributions

# Input format can be either averaged or replicates listed side-by-side; units are arbitrary (counts, fpkm, vsd, rlog, etc) 

# SNF Winter 2014/2015

import sys, os, argparse, csv, random, re
from SNFUtils import * # prompt, safe_open_file, stdout_from_command
from collections import defaultdict, namedtuple
from scipy.stats.stats import pearsonr, spearmanr, kurtosis
from scipy.optimize import curve_fit
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt 
from itertools import izip_longest

parser = argparse.ArgumentParser(description="Plot input transcript ID from input distribution file") 

#parser.add_argument("-i", "--input", help="File containing transcript expression values", required=True)
parser.add_argument("-i", "--input", help="File containing transcript distributions", required=True)
parser.add_argument("-o", "--output", help="Output filename (default is stdout)", default="stdout")
parser.add_argument("-n", "--nrep", help="Number of replicates of each point", required=True, type=int)
parser.add_argument("--id", help="Transcript ID(s) to print (can be partial; can be comma-separated list)", required=True) 
parser.add_argument("--tx-to-gene", help="File containing transcript ID to gene name mapping", required=True) 
parser.add_argument("--text", help="Output text data in addition to plots.", action="store_true") 

args = parser.parse_args()

print "----------------------------------------"
print "|     plot_tripseq_transcript.py        |"
print "|   plot transcript distribution        |"
print "|        run with -h for help           |"
print "|             snf 2/2015                |"
print "------------------------------------\n\n"

print "Arguments:"
print args

if (args.output == "stdout"):
    outfile = sys.stdout
else:
    outfile = safe_open_file(args.output)

tx_to_gene = defaultdict(str)
tx_to_name = defaultdict(str)
tx_to_txname = defaultdict(str)

processed = 0

print "\nCreating tx_to_gene dictionary from %s..." % args.tx_to_gene

with open(args.tx_to_gene, "r") as txfile:
    for line in txfile:
        line = line.split()
        if (len(line) < 4): 
            sys.exit("ERROR: line in tx_to_gene file has less than four entries: %s" % line)

        tx_to_gene[line[0]] = line[1]
        tx_to_name[line[0]] = line[2]
        tx_to_txname[line[0]] = line[3]

        processed += 1

print "Dictionary created (%d entries)." % processed

# this is used to define the order of the input columns in the output plots, to rearrange input to output.
# janky but it works for these input files. 
ORDER = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 3]

ids = args.id.split(",") 

ids = [ ids[i].strip() for i in range(len(ids))]

print "Looking for ids %s" % ids 

# sample input format:

# -- first column must be transcript ID; quotes are optional
# -- header line below is optional and detected by a non-numeric second column
# "ensid" "X80S"  "poly2" "poly3" "poly4" "poly5" "poly6" "poly7" "poly8"
#"ENST00000000233"       -0.562625014929129      -2.23696624248826       -0.765273156804401      0.702668471793972       1.36603995862049        1.02457718119074        0.751424784189265       -0.279845981572665
#"ENST00000000412"       -1.10133854138933       -1.77019107516088       -1.7623840822903        0.404249332913166       1.42491481379874        1.42016301836686        1.35747792898968        0.0271086047720868
#"ENST00000000442"       1.1367571728999 -2.58192112517174       0.0240875688746112      0.482589702447365       0.969885993274519       0.643070451106668       -0.177572605307724      -0.496897158123603

np.seterr(all='raise')

processed = 0
rejected_pe = 0
rejected_zero = 0

fname = args.input

with open(fname, "r") as infile:
    #inspect the first line
        
    random.seed(1)
        
    if fname.endswith("csv"):
        line = [i.strip("\"") for i in infile.readline().split(",")]
    else: 
        line = [i.strip("\"") for i in infile.readline().split()]

    if (is_number(line[2])): # no header
        found_header = False
        infile.seek(0) #rewind
        print "\tNo header found."

    else:
        found_header = True
        header = line
        print "\tHeader found on line one."
        print "\tOutput order (CHECK THIS; should be 40, 60, p1-8, cyto):"
        header = [ header[i*args.nrep+1] for i in ORDER ]
        print "\t\t" + ", ".join(header)

    if ( (len(line) - 1) % args.nrep):
        print "\tWARNING: %d observations found; not divisible by %d replicates." % (len(line) - 1, args.nrep) 


    npts = (len(line) - 1)/args.nrep
    print "\tProcessing %d points (%d replicates)" % ( npts, args.nrep)

                                                           


    for line in infile:
            
        if (not line.strip()):
            continue

        if fname.endswith("csv"):
            line = [i.strip("\"") for i in line.split(",")]
        else: 
            line = [i.strip("\"") for i in line.split()]

            # only process the one transcript 
        if True in [ re.search(ids[i], line[0]) != None for i in range(len(ids))] or True in [tx_to_name[line[0]] == ids[i] for i in range(len(ids))]:
            print "\tPlotting id %s" % line[0]

                # convert all to floats before continuing...
            line = [line[0]] + [float(i) for i in line[1:]]

            txid = line[0]
            obs = line[1:]  # the observations
            nobs = len(obs) / args.nrep
  
                # 1 read on a 1kb transcript out of 10 million -> FPKM = 1e-10.  rounding to this precision.  
            obs = np.around(obs, 10)
            
            if (len(obs) % args.nrep): 
                print "\tWarning: transcript %s has %d observations which is not divisible by %d replicates" % (txid, len(obs), args.nrep)

                # collapse replicates into means, preserve stdev
            means = [np.mean(i) for i in chunkwise(obs, args.nrep)]
            stdevs = [np.std(i) for i in chunkwise(obs, args.nrep)]

            # reorder the lists here to go 40S, 60S, 80S, p2-p8, cyto
            means = [ means[i] for i in ORDER ]
            stdevs = [ stdevs[i] for i in ORDER ] 

            if txid in tx_to_name:
                #plotname = txid + " (" + tx_to_name[txid] + ")"
                plotname = tx_to_txname[txid] + " (" + txid + ")"
                outfilename = tx_to_txname[txid] + "_" + txid + ".pdf"
                if (args.text):
                    outdataname = tx_to_txname[txid] + "_" + txid + ".txt" 
            else:
                print "Warning: gene name for transcript %s not found" % txid
                plotname = txid
                outfilename = plotname + ".pdf"
                if (args.text):
                    outdataname = plotname + ".txt" 

            plot_dist_fancy(range(1,npts+1), means, outfilename, stdevs, plotname)

            if (args.text):
                outfile = safe_open_file(outdataname)
            
                outfile.write("\n".join( [ "%s\t%d\t%.6f\t%.6f" % (header[i], i, means[i], stdevs[i]) for i in range(len(header))]))
                outfile.write("\n")

