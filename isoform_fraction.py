#!/usr/bin/env python

# isoform_fraction.py

# calculate the percent of genes with a "minor" isoform above a specified fraction or range of fractions

# SNF 2017

import sys, os, argparse, csv, random, re
from SNFUtils import * # prompt, safe_open_file, stdout_from_command
from collections import defaultdict, namedtuple
from scipy.stats import mannwhitneyu, ks_2samp
import numpy as np
import matplotlib.pyplot as plt
from itertools import izip_longest

parser = argparse.ArgumentParser(description="Calculate the percent of genes with a minor isoform above a specified fraction")

parser.add_argument("-i", "--input", help="File containing transcript distributions", required=True)
parser.add_argument("-o", "--output", help="Output filename (default is stdout)", default="stdout")
parser.add_argument("-n", "--nrep", help="Number of replicates of each point", required=True, type=int)
parser.add_argument("--tx-to-gene", help="File containing transcript ID to gene name mapping", required=True)
parser.add_argument("--text", help="Output text data in addition to plots.", action="store_true")
parser.add_argument("--format", help="Image format to export (png or pdf).", default="pdf")
parser.add_argument("--isoform-fraction", help="What fraction of minor isoform to include (default 0.1)", type=float, default=0.1)


args = parser.parse_args()

if (args.output == "stdout"):
    outfile = sys.stdout
else:
    outfile = safe_open_file(args.output)

tx_to_gene = defaultdict(str)
tx_to_name = defaultdict(str)
tx_to_txname = defaultdict(str)

print "Assuming five samples: cyto, nuc, monosome, low polysome, high polysome."
Isoform_samples = namedtuple('Isoform_samples', ['id', 'cyto', 'nuc', 'mono', 'low', 'high'])

# each of these is a dictionary with keys of gene names with values a list namedtuple of Samples for each isoform
hesc = defaultdict(list)
npc = defaultdict(list)
neu14 = defaultdict(list) #14 day neurons
neu50 = defaultdict(list) #50 day neurons

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

# sample input format:


#txid s1r1 s1r2 s2r1 s2r2 s3r1 s3r2 etc

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

        # cufflinks output has the replicate number as the last two (i.e. _0 or _1). this trims.
        header = [i[:-2] for i in line]

        print "\tHeader found on line one."
        print "\tOrder parsed: "

        header = [ header[i*args.nrep+1] for i in range(len(header)/2) ]

        print "\t\t" + ", ".join(header)
        print "Header has %d columns" % len(header)

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
        #stdevs = [np.std(i) for i in chunkwise(obs, args.nrep)]

        if txid in tx_to_gene:
            geneid = tx_to_gene[txid]
            hesc[geneid].append(Isoform_samples(txid, *means[0:5]))
            npc[geneid].append(Isoform_samples(txid, *means[5:10]))
            neu14[geneid].append(Isoform_samples(txid, *means[10:15]))
            neu50[geneid].append(Isoform_samples(txid, *means[15:20]))

        else:
            print "Warning: gene name for transcript %s not found" % txid

        processed += 1

        if not (processed % 5000):
            print "Processed %d lines." % processed

# now, have read in all the transcripts into these gene dicts. compute the percent that have above whatever fraction

# calculate the percent of genes that have a minor transcript above the specified fraction for dicts constructed as described above
# returns this percent as a fraction

def calc_isoform_fraction(dict, isoform_fraction, sample):
    genes_over_fraction = 0
    genes_processed = 0
    geneids_over_fraction = []
    samples = []

    for key, val in dict.iteritems():
        if sample == "cyto":
            samples = [i.cyto for i in val]
        elif sample == "nuc":
            samples = [i.nuc for i in val]
        elif sample == "mono":
            samples = [i.mono for i in val]
        elif sample == "low":
            samples = [i.low for i in val]
        elif sample == "high":
            samples = [i.high for i in val]
        else:
            print "Unrecognized sample %s" % sample

        total = sum(samples)
        if total > 1 : # only continue if expressed (TPM > 1 for the gen) 

        # transform to fractional isoform abundance
            samples = [i / total for i in samples]

            n_over_fraction = sum([i > isoform_fraction for i in samples])

            if n_over_fraction > 1: # implies at least one minor isoform above the fraction cutoff
                genes_over_fraction += 1
                geneids_over_fraction.append(key)

            genes_processed += 1

    #print "Found %d genes over fraction %.2f out of %d genes (%3.2f%%) in %s" % (genes_over_fraction, isoform_fraction, genes_processed, genes_over_fraction/float(genes_processed)*100, sample)

    #print geneids_over_fraction[0:10]


    return genes_over_fraction/float(genes_processed)

#for sample in ["cyto", "nuc", "mono", "low", "high"]:
#    calc_isoform_fraction(hesc, args.isoform_fraction, sample)
#    calc_isoform_fraction(npc, args.isoform_fraction, sample)
#    calc_isoform_fraction(neu14, args.isoform_fraction, sample)
#    calc_isoform_fraction(neu50, args.isoform_fraction, sample)

def plot_it(cell_dict, name, npoints): 

    xvals = np.linspace(0, 0.5, npoints)
    cyto_y = [calc_isoform_fraction(cell_dict, xvals[i], "cyto") for i in range(len(xvals))]
    nuc_y = [calc_isoform_fraction(cell_dict, xvals[i], "nuc") for i in range(len(xvals))]
    mono_y = [calc_isoform_fraction(cell_dict, xvals[i], "mono") for i in range(len(xvals))]
    low_y = [calc_isoform_fraction(cell_dict, xvals[i], "low") for i in range(len(xvals))]
    high_y = [calc_isoform_fraction(cell_dict, xvals[i], "high") for i in range(len(xvals))]

    cdffile = safe_open_file("%s_ecdf.csv" % name)



    cdffile.write("xvals,%s_cyto,%s_nuc,%s_mono,%s_low,%s_high\n" % (name, name, name, name, name))
    cdffile.write("\n".join(["%f,%f,%f,%f,%f,%f" % (xvals[i], cyto_y[i], nuc_y[i], mono_y[i], low_y[i], high_y[i]) for i in range(len(xvals))]))
    
    mw_cn_test = mannwhitneyu(cyto_y, nuc_y)
    mw_ch_test = mannwhitneyu(cyto_y, high_y)

    # plot eCDFs
    fig, ax = plt.subplots()
    
    plt.plot( xvals, cyto_y )
    plt.plot( xvals, nuc_y )
    plt.plot( xvals, mono_y )
    plt.plot( xvals, low_y )
    plt.plot( xvals, high_y )

    plt.legend(["cyto", "nuc", "mono", "low", "high"], fontsize=8)

    #plt.text(.7, 0.15, "K-S, p = %4.2f, %.3g" % (ks_test[0], ks_test[1]), fontsize=8, transform=ax.transAxes)
    plt.text(.7, 0.15, "M-W CN U, p = %4.2f, %.3g" % (mw_cn_test[0], mw_cn_test[1]), fontsize=8, transform=ax.transAxes)
    plt.text(.7, 0.1, "M-W CH U, p = %4.2f, %.3g" % (mw_ch_test[0], mw_ch_test[1]), fontsize=8, transform=ax.transAxes)
#plt.text(.7, 0.05, "set1: %s; set2: %s" % (", ".join(args.set1[1::2]), ", ".join(args.set2[1::2])), fontsize=8, transform=ax.transAxes)
#plt.text(.7, 0.02, "n_1: %d; n_2: %d; n_3: %d" % (len(ivals), len(jvals)), fontsize=8, transform=ax.transAxes)

    plt.savefig("%s_ecdf.png" % name)
    plt.close(fig)

print "plotting hesc"
plot_it(hesc, "hesc", 100)
print "plotting npc"
plot_it(npc, "npc", 100) 
print "plotting neu14"
plot_it(neu14, "neu14", 100)
print "plotting neu50"
plot_it(neu50, "neu50", 100) 


# now count how many isoforms there are per gene that are above a certain fraction 

# returns a dict of [number_of_isoforms][count_of_genes]

def isoforms_per_gene(dict, fraction, sample): 
    isoforms_over_fraction = defaultdict(int) 
    
    samples = []

    for key, val in dict.iteritems():
        if sample == "cyto":
            samples = [i.cyto for i in val]
        elif sample == "nuc":
            samples = [i.nuc for i in val]
        elif sample == "mono":
            samples = [i.mono for i in val]
        elif sample == "low":
            samples = [i.low for i in val]
        elif sample == "high":
            samples = [i.high for i in val]
        else:
            print "Unrecognized sample %s" % sample

        total = sum(samples)
        if total > 1: # only continue if expressed (TPM > 1 cutoff)

        # transform to fractional isoform abundance
            samples = [i / total for i in samples]

            n_over_fraction = sum([i > fraction for i in samples])
            
            isoforms_over_fraction[n_over_fraction] += 1

    return isoforms_over_fraction

# run isoforms_per_gene and print the output for a cell type (name) for a specified percent (not fraction)

def count_it(dict, name, percent): 

    cyto = isoforms_per_gene(dict, percent/100., "cyto") 
    nuc = isoforms_per_gene(dict, percent/100., "nuc") 
    mono = isoforms_per_gene(dict, percent/100., "mono") 
    low = isoforms_per_gene(dict, percent/100., "low") 
    high = isoforms_per_gene(dict, percent/100., "high") 

    histfile = safe_open_file("%s_isoforms_%d.csv" % (name, percent))

    nrow = max(max(cyto.keys()), max(nuc.keys()), max(mono.keys()), max(low.keys()), max(high.keys()))
    xvals = range(nrow+1)

    cyto_list = [cyto[i] for i in xvals]
    nuc_list = [nuc[i] for i in xvals]
    mono_list = [mono[i] for i in xvals]
    low_list = [low[i] for i in xvals]
    high_list = [high[i] for i in xvals]

    histfile.write("n_isoforms,%s_cyto,%s_nuc,%s_mono,%s_low,%s_high\n" % (name, name, name, name, name))
    histfile.write("\n".join(["%f,%f,%f,%f,%f,%f" % (xvals[i], cyto[i], nuc[i], mono[i], low[i], high[i]) for i in range(nrow+1)]))

    # plot "histograms" 
    fig, ax = plt.subplots()
    

    plt.plot( xvals, cyto_list )
    plt.plot( xvals, nuc_list )
    plt.plot( xvals, mono_list )
    plt.plot( xvals, low_list )
    plt.plot( xvals, high_list )

    plt.legend(["cyto", "nuc", "mono", "low", "high"], fontsize=8)

    plt.savefig("%s_hist_%d.png" % (name, percent))
    plt.close(fig)


for i in [5, 10, 25]:
    print "counting hesc %d" % i
    count_it(hesc, "hesc", i)
    print "counting npc %d" % i
    count_it(npc, "npc", i)
    print "counting neu14 %d" % i
    count_it(neu14, "neu14", i)
    print "counting neu50 %d" % i
    count_it(neu50, "neu50", i)
