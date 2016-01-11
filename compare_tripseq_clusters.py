#!/usr/bin/env python

# compare_tripseq_clusters.py

# Compare input lists of clusters against each other searching for genes that have transcripts in the different clusters.

# Intended to be applied to fractionated polysome sequencing data (a la tripseq) but theoretically can work for any data involving transcript-specific distributions

# Input format can be either averaged or replicates listed side-by-side; units are arbitrary (counts, fpkm, vsd, rlog, etc)

# SNF Winter 2014/2015

import sys, os, argparse, csv, random
import GTF
from Transcript import *
from SNFUtils import * # prompt, safe_open_file, stdout_from_command
from collections import defaultdict, namedtuple
from scipy.stats.stats import pearsonr, spearmanr, kurtosis
from scipy.optimize import curve_fit
from scipy.stats import mannwhitneyu, ks_2samp
import numpy as np
import matplotlib.pyplot as plt
from itertools import izip_longest

def poly3(x, p1, p2, p3, p4):
    return p1 + p2 * x + p3 * x**2 + p4 * x**3

parser = argparse.ArgumentParser(description="Compare transcript-specific distributions at the gene level")

#parser.add_argument("-i", "--input", help="File containing transcript expression values", required=True)
parser.add_argument("--set1", help="Files and IDs containing transcript distributions; compare between set1 and set2", metavar="FNAME ID ...", nargs="+", required=True)
parser.add_argument("--set2", help="Files and IDs containing transcript distributions; compare between set1 and set2", metavar="FNAME ID ...", nargs="+", required=True)
parser.add_argument("--tx-to-gene", help="Mapping between transcript ID in input file and gene ID", required=True)
parser.add_argument("-o", "--output", help="Output filename (default is stdout)", default="stdout")
parser.add_argument("-n", "--nrep", help="Number of replicates of each point", required=True, type=int)
parser.add_argument("--txome-props", help="List of files with transcriptome properties to correlate among (wildcards ok)", nargs="+", default=None)
parser.add_argument("--control", help="Perform randomized comparisons of input transcripts as a control.", action="store_true")
parser.add_argument("--txome-gtf", help="Path to transcriptome GTF", required=True)
parser.add_argument("--strip-id", help="Strip anything following the final _ off the transcript ID in input", action="store_true") 

# below not implemented yet
#parser.add_argument("--seq", help="Output sequence of the region(s) different between transcripts", action="store_true")
#parser.add_argument("-g", "--genome", help="Path to the genome file (required for --seq")

args = parser.parse_args()

print "------------------------------------"
print "|    compare_tripseq_clusters.py   |"
print "| extract differential transcripts |"
print "|      run with -h for help        |"
print "|           snf 1/2015             |"
print "------------------------------------\n\n"

print "Arguments:"
print args

# sanity check inputs...
if (len(args.set1) % 2 or len(args.set2) % 2):
    sys.exit("FATAL: --set1 and --set2 must be lists of file id file id ...")

if (args.control):
    print "\nPerforming a comparsion between two randomized subsets of the data as a control"

if (args.output == "stdout"):
    outfile = sys.stdout
else:
    outfile = safe_open_file(args.output)

# First, create a dictionary mapping from transcript ID to gene ID by reading in args.tx_to_gene
#  - optionally, also save the transcript's gene name and biotype if they are present in columns 3 and 4 in the input
tx_to_gene = defaultdict(str)
tx_to_name = defaultdict(str)
tx_to_biotype = defaultdict(str)

processed = 0

print "\nCreating tx_to_gene dictionary from %s..." % args.tx_to_gene

with open(args.tx_to_gene, "r") as txfile:
    for line in txfile:
        line = line.split()
        tx_to_gene[line[0]] = line[1]

        if (len(line) > 2):
            tx_to_name[line[0]] = line[2]

        if (len(line) > 3):
            tx_to_biotype[line[0]] = line[3]

        processed += 1

print "Dictionary created (%d entries)." % processed

print "\nCreating transcript dictionary from %s (currently commented out)" % args.txome_gtf

txDict = defaultdict(list)
processed = 0

# COMMENTED CODE builds a transcript dictionary to extract transcripts with specific differences (i.e. 5' utr changes, etc).  this takes a bit, so commented for now.

# the issue here is that lines for various transcripts may be interleaved, so can either create lots of objects, or a giant dict. opted for giant dict.
# for line in GTF.lines(args.txome_gtf):

#     txDict[line["transcript_id"]].append(line)
#     processed += 1

#     if (not processed % 25000):
#         print "\tProcessed %d lines..." %  processed

# print "Dictionary built (%d entries)." % processed

# print "Creating transcript objects."
# processed = 0

# # now create a Transcript object for each transcript and output it

# for key in txDict:

#     tx = createGTFTranscript(txDict[key])

#     txDict[key] = tx
#         #---------------- save to a new dict and/or overwrite old one here

#             #print tx
#         #writeOutput(tx)
#     processed += 1

#     if (not processed % 2500):
#         print "\tProcessed %d entries..." %  processed

print "Dictionary built (%d entries)." % processed
# sample input format:

# -- first column must be transcript ID; quotes are optional
# -- header line below is optional and detected by a non-numeric second column
# "ensid" "X80S"  "poly2" "poly3" "poly4" "poly5" "poly6" "poly7" "poly8"
#"ENST00000000233"       -0.562625014929129      -2.23696624248826       -0.765273156804401      0.702668471793972       1.36603995862049        1.02457718119074        0.751424784189265       -0.279845981572665
#"ENST00000000412"       -1.10133854138933       -1.77019107516088       -1.7623840822903        0.404249332913166       1.42491481379874        1.42016301836686        1.35747792898968        0.0271086047720868
#"ENST00000000442"       1.1367571728999 -2.58192112517174       0.0240875688746112      0.482589702447365       0.969885993274519       0.643070451106668       -0.177572605307724      -0.496897158123603

np.seterr(all='raise')

print "\nProcessing input transcripts..."

def process_file(fname, id, gene_to_txs, tx_to_fit):
    processed = 0
    rejected_pe = 0
    rejected_zero = 0

    with open(fname, "r") as infile:
    #inspect the first line

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

        if ( (len(line) - 1) % args.nrep):
            print "\tWARNING: %d observations found; not divisible by %d replicates." % (len(line) - 1, args.nrep)

        print "\tProcessing %d points (%d replicates)" % ( (len(line) - 1)/args.nrep, args.nrep)




        for line in infile:

            if (not line.strip()):
                continue

            if fname.endswith("csv"):
                line = [i.strip("\"") for i in line.split(",")]
            else:
                line = [i.strip("\"") for i in line.split()]

            # append all transcripts for a particular gene to a dictionary of lists
            if (tx_to_gene[line[0]]):

                # convert all to floats before continuing...
                line = [line[0]] + [float(i) for i in line[1:]]

                txid = line[0]
                geneid = tx_to_gene[txid]
                obs = line[1:]  # the observations
                nobs = len(obs) / args.nrep

                # 1 read on a 1kb transcript out of 10 million -> FPKM = 1e-10.  rounding to this precision.
                obs = np.around(obs, 10)

                if (len(obs) % args.nrep):
                    print "\tWarning: transcript %s has %d observations which is not divisible by %d replicates" % (txid, len(obs), args.nrep)

                # collapse replicates into means, preserve stdev
                means = [np.mean(i) for i in chunkwise(obs, args.nrep)]
                stdevs = [np.std(i) for i in chunkwise(obs, args.nrep)]

                # average percent error

                mean_pe = np.mean([stdevs[i]/means[i] for i in range(nobs) if means[i] > 1e-6])

                # if more than half are zero; reject
                if ( sum( [np.abs(i) < 0.01 for i in means] ) > nobs/2):
                    #print "Disqualifying %s: more than half the values are zero (%s)" % (txid, means)
                #third_order_poly_fit_plot(range(1,nobs+1), means, txid + "_stdevreject.png", stdevs)
                    rejected_zero += 1


            # if this tx wasn't disqualified, add it to the mix
                else:
                    gene_to_txs[geneid].append( ([txid] + list(obs), [txid] + means, [txid] + stdevs, id) )
                    # save the fit for later.
                    #tx_to_fit[txid] = (popt, pcov, infodict, mesg, ier, xvals, means, xnew, ynew)


                processed += 1
                if (not processed % 5000):
                    print "\tProcessed %d entries..." % processed
            # plot the transcript..
            # nasty list comprehensions to make a list of nrep repetitions of each obs to generate the x coordinate list
            #xvals = [ k for j in [[i] * args.nrep for i in range(1,nobs+1)] for k in j]
            #third_order_poly_fit_plot(xvals, obs, txid + "_poly.png", None)
            #third_order_poly_fit_plot(range(1,nobs+1), means, txid + "_poly_avg.png", stdevs)

            else:
                print "\tWarning: gene ID not found for transcript ID %s; skipping." % line[0]
        #print line

    #total = float(processed + rejected_pe + rejected_zero )

    print "\tRead %d transcript distributions from %s; rejected %d (%.1f %%) for percent error; %d (%.1f %%) for half zero (total %.1f %%)" \
        % (processed, fname, rejected_pe, rejected_pe/float(processed)*100, rejected_zero, rejected_zero/float(processed)*100, (rejected_pe + rejected_zero)/float(processed)*100)
    return processed - rejected_pe - rejected_zero
    #return (gene_to_txs, tx_to_fit)


# workhorse data structures for the program - holds all transcript distributions as values for the gene (key)
gene_to_txs_set1 = defaultdict(list)
tx_to_fit_set1 = defaultdict(tuple)
gene_to_txs_set2 = defaultdict(list)
tx_to_fit_set2 = defaultdict(tuple)

#gene_to_txs_mean = defaultdict(list)
#gene_to_txs_stdev = defaultdict(list)

# global indices into the dict of lists
IDXDIST = 0
IDXMEAN = 1
IDXSTDEV = 2
IDXCLUSTID = 3

n_in_set1 = 0

for fname,id in pairwise(args.set1):
    print "\nReading from set1 file %s [ID %s]" % (fname, id)
    n_in_set1 += process_file(fname, id, gene_to_txs_set1, tx_to_fit_set1)

print "Read %d total transcripts for set1." % n_in_set1

n_in_set2 = 0

for fname,id in pairwise(args.set2):
    print "\nReading from set2 file %s [ID %s]" % (fname, id)
    n_in_set2 += process_file(fname, id, gene_to_txs_set2, tx_to_fit_set2)

print "Read %d total transcripts for set2." % n_in_set2

print "%d genes in set1, %d genes in set2" % (len(gene_to_txs_set1), len(gene_to_txs_set2))

#print "---- gene_to_txs_set1 ----"
#print gene_to_txs_set1
#print tx_to_fit_set1

# read in transcriptome parameters:

txome_props = {}

#with open("/Volumes/11tb/seq/Homo_sapiens/Ensembl/GRCh37/Annotation/snf_transcriptome_annotations/grch37_cds_deltag_cds.csv", "r") as infile:
#    reader = csv.reader(infile)
#    mydict = {rows[0]:rows[1] for rows in reader}

def row_to_float(row):
    return [float(i.strip("\";")) if is_number(i.strip("\";")) else i.strip("\";") for i in row[1:]]
    #if is_number(row[1]):
    #    return np.array(row[1:], dtype=np.float32)
    #else:
    #    return row[1:]

if (args.txome_props):
    if "*" in args.txome_props:
        files = glob.glob(args.txome_props)
    else:
        files = args.txome_props

    for file in files:
        print "Reading transcriptome properties from %s" % file
        with open(file, "r") as infile:
            reader = csv.reader(infile)
            header = reader.next()
            prop_name = os.path.splitext(os.path.basename(file))[0]
            if (args.strip_id): 
                txome_props[prop_name] = {"_".join(rows[0].split("_")[:-1]).strip("\";"):row_to_float(rows)  for rows in reader}
            else: 
                txome_props[prop_name] = {rows[0].strip("\";"):row_to_float(rows)  for rows in reader}
                
            #txome_props[prop_name] = {rows[0].split('_')[0]:row_to_float(rows)  for rows in reader}
            txome_props[prop_name]["header"] = header[1:]
            print "\tRead %d entries." % len(txome_props[prop_name].keys())

# generate two randomized, size-matched sets as a control, if requested

gene_to_txs_rand1 = defaultdict(list)
gene_to_txs_rand2 = defaultdict(list)
n_in_rand1, n_in_rand2 = 0,0

if (args.control):
    set1_txs = gene_to_txs_set1.viewvalues()
    set2_txs = gene_to_txs_set2.viewvalues()

    # values are list of tuples of (dist, mean, stdev, clusterID) where each list item is a transcript.
    # indices are IDXDIST - distribution, IDXMEAN - means, IDXSTDEV - stdevs, IDXCLUSTID - cluster identifier

    # collapse lists and combine...
    all_txs = [k for j in list(set1_txs) for k in j] + [k for j in list(set2_txs) for k in j]
    print "\nGenerating random transcripts from list of %d (set1 has %d, set2 has %d)" % (len(all_txs), n_in_set1, n_in_set2)

    random.seed(1)
    all_txs = random.sample(all_txs, len(all_txs))

    for tx in all_txs[:n_in_set1]:
        txid = tx[IDXDIST][0]
        gene_to_txs_rand1[tx_to_gene[txid]].append(tx)
        n_in_rand1 += 1

    for tx in all_txs[n_in_set1:]:
        txid = tx[IDXDIST][0]
        gene_to_txs_rand2[tx_to_gene[txid]].append(tx)
        n_in_rand2 += 1


    #print "set1: %d rand1: %d set2: %d rand2: %d" % (n_in_set1, n_in_rand1, n_in_set2, n_in_rand2)

# now, for all genes that have transcripts in set1 and set2, report things.

print "\nComparing transcript parameters."

prop_list_set1 = defaultdict(list)
prop_list_set2 = defaultdict(list)
tx_list_set1 = defaultdict(list)
tx_list_set2 = defaultdict(list)

comparisons, length_different = 0,0
length_ratios = []

genes_in_both_sets, txs_compared_in_set1, txs_compared_in_set2 = 0,0,0

for gene in gene_to_txs_set1.iterkeys():

    # this data structure gene_to_txs is now for each gene (gene), a list (txs) of tuples of (dist, mean, stdev, clusterID) where each entry in the list is one transcript.
    # possible that genes will have multiple transcripts in both set1 and set2 - need to do all-by-all comparisons when reporting things

    # indices into the tuples are IDXDIST - distribution, IDXMEAN - means, IDXSTDEV - stdevs, IDXCLUSTID - cluster identifier

    if gene in gene_to_txs_set2:
        genes_in_both_sets += 1
        # catalog properties for the two sets
        for val1 in gene_to_txs_set1[gene]:
            txs_compared_in_set1 += 1
            tx_i = val1[IDXDIST][0]
            for propfile, propdict in txome_props.iteritems():
                for k in range(len(propdict["header"])):
                    propname = propfile + "_" + propdict["header"][k]
                    if tx_i in propdict:
                        prop_i = propdict[tx_i][k]
                        prop_list_set1[propname].append(prop_i)
# current construction of this does not handle missing values
                        tx_list_set1[tx_i].append(str(prop_i))
                    else:
                        tx_list_set1[tx_i].append(str(-1))
                    #else:
                    #    print "Warning: transcript %s not found for property %s" % (tx_i, propfile)

        for val2 in gene_to_txs_set2[gene]:
            txs_compared_in_set2 += 1
            tx_j = val2[IDXDIST][0]
            for propfile, propdict in txome_props.iteritems():
                for k in range(len(propdict["header"])):
                    propname = propfile + "_" + propdict["header"][k]
                    if tx_j in propdict:
                        prop_j = propdict[tx_j][k]
                        prop_list_set2[propname].append(prop_j)
                        tx_list_set2[tx_j].append(str(prop_j))
                    else:
                        tx_list_set2[tx_j].append(str(-1))
                    #else:
                    #    print "Warning: transcript %s not found for property %s" % (tx_j, propfile)

        # all by all length comparison...

        # this is hackish and assumes certain input property filenames AND that set1 is high polys and set2 is low polys

        for val1 in gene_to_txs_set1[gene]:
            for val2 in gene_to_txs_set2[gene]:
                tx1 = val1[IDXDIST][0]
                tx2 = val2[IDXDIST][0]

                comparisons += 1

                # propdict = txome_props["grch37_cds_length"]
#                 if tx1 in propdict and tx2 in propdict:
#                     len1 = propdict[tx1][0]
#                     len2 = propdict[tx2][0]

#                     length_ratios.append(len2/float(len1))

#                     if (len2 > len1):
#                         #print "For gene %s tx1 %s length (%d) is < tx2 %s length (%d)!" % (gene, tx1, len1, tx2, len2)
#                         length_different += 1

#print "Compared %d transcripts, %d have longer set2 length than set1 length (%3.2f %%)" % (comparisons, length_different, length_different/float(comparisons)*100)
print "Compared %d genes" % genes_in_both_sets
print "Compared %d transcripts overall, %d from set1 and %d from set2" % (comparisons, txs_compared_in_set1, txs_compared_in_set2)

# plot histogram of the length ratios
# fig, ax = plt.subplots()
# n, bins, patches = plt.hist(length_ratios, bins=np.logspace(-1,1,100), range=(0.1,10), normed=1, facecolor='green', alpha=0.75)
# ax.set_xscale("log")
# plt.xlabel("len_set2/len_set1")
# plt.savefig("length_ratio_hist.png")
# plt.close(fig)


# Generate random controls -------

prop_list_rand1 = defaultdict(list)
prop_list_rand2 = defaultdict(list)

comparisons, length_different = 0,0

if (args.control):
    tested = 0

    for gene in gene_to_txs_rand1.iterkeys():

    # this data structure gene_to_txs is now for each gene (gene), a list (txs) of tuples of (dist, mean, stdev, clusterID) where each entry in the list is one transcript.
    # possible that genes will have multiple transcripts in both rand1 and rand2 - need to do all-by-all comparisons when reporting things

    # indices into the tuples are IDXDIST - distribution, IDXMEAN - means, IDXSTDEV - stdevs, IDXCLUSTID - cluster identifier

        if gene in gene_to_txs_rand2:
            # catalog properties for the two rands
            for val1 in gene_to_txs_rand1[gene]:
                tx_i = val1[IDXDIST][0]
                for propfile, propdict in txome_props.iteritems():
                    for k in range(len(propdict["header"])):
                        propname = propfile + "_" + propdict["header"][k]
                        if tx_i in propdict:
                            prop_i = propdict[tx_i][k]
                            prop_list_rand1[propname].append(prop_i)
                        #else:
                        #    print "Warning: transcript %s not found for property %s" % (tx_i, propfile)

            for val2 in gene_to_txs_rand2[gene]:
                tx_j = val2[IDXDIST][0]
                for propfile, propdict in txome_props.iteritems():
                    for k in range(len(propdict["header"])):
                        propname = propfile + "_" + propdict["header"][k]
                        if tx_j in propdict:
                            prop_j = propdict[tx_j][k]
                            prop_list_rand2[propname].append(prop_j)
                        #else:
                        #    print "Warning: transcript %s not found for property %s" % (tx_j, propfile)
            for val1 in gene_to_txs_set1[gene]:
                for val2 in gene_to_txs_set2[gene]:
                    tx1 = val1[IDXDIST][0]
                    tx2 = val2[IDXDIST][0]

                    comparisons += 1

#                     propdict = txome_props["grch37_cds_length"]
#                     if tx1 in propdict and tx2 in propdict:
#                         len1 = propdict[tx1][0]
#                         len2 = propdict[tx2][0]

#                         length_ratios.append(len2/float(len1))

#                         if (len2 > len1):
#                         #print "For gene %s tx1 %s length (%d) is < tx2 %s length (%d)!" % (gene, tx1, len1, tx2, len2)
#                             length_different += 1

#print "Compared %d randomly assigned transcripts, %d have longer set2 length than set1 length (%3.2f %%)" % (comparisons, length_different, length_different/float(comparisons)*100)
print "Compared %d randomly assigned transcripts" % (comparisons)

# now make eCDF plots...

def mkline(tuple):
    if False in [is_number(i) for i in tuple]:
        return ",".join(["%s" % i if i is not None else "" for i in tuple])
    elif str in [type(i) for i in tuple]:
        print "STR DETECT %s " % str(tuple)
    else:
        return ",".join(["%.6f" % i if i is not None else "" for i in tuple])

for propname in prop_list_set1.iterkeys():

    print "Plotting property %s" % propname

    ivals = prop_list_set1[propname]
    jvals = prop_list_set2[propname]

    propfile = safe_open_file("%s_ecdf.csv" % propname)

    if args.control:
        rivals = prop_list_rand1[propname]
        rjvals = prop_list_rand2[propname]

    #if (len(ivals) != len(jvals)):
    #    print "Warning: property %s i and j lists are unequal length." % propname

    if False in [is_number(i) for i in ivals] or False in [is_number(j) for j in jvals]:
        print "Skipping property %s due to non-numeric values" % propname

        if args.control:
            propfile.write("set1,set2,rand1,rand2\n")
            propfile.write("\n".join([mkline(i) for i in izip_longest(ivals, jvals, rivals, rjvals)]))
            propfile.close()

        else:
            propfile.write("set1,set2\n")
            propfile.write("\n".join([mkline(i) for i in izip_longest(ivals, jvals)]))
            propfile.close()

        continue

    # calculate Mann Whitney U p-value
    mw_test = mannwhitneyu(ivals, jvals)

    # calculate the KS statistics
    ks_test = ks_2samp(ivals, jvals)

    # plot eCDFs
    fig, ax = plt.subplots()

    i_sorted=np.sort( ivals )
    i_yvals=np.arange(len(i_sorted))/float(len(i_sorted))
    plt.plot( i_sorted, i_yvals )

    j_sorted=np.sort( jvals )
    j_yvals=np.arange(len(j_sorted))/float(len(j_sorted))
    plt.plot( j_sorted, j_yvals )
    plt.legend([propname + "_set1", propname + "_set2"], fontsize=8)

    plt.text(.7, 0.15, "K-S, p = %4.2f, %.3g" % (ks_test[0], ks_test[1]), fontsize=8, transform=ax.transAxes)
    plt.text(.7, 0.1, "M-W U, p = %4.2f, %.3g" % (mw_test[0], mw_test[1]), fontsize=8, transform=ax.transAxes)
    plt.text(.7, 0.05, "set1: %s; set2: %s" % (", ".join(args.set1[1::2]), ", ".join(args.set2[1::2])), fontsize=8, transform=ax.transAxes)
    plt.text(.7, 0.02, "n_1: %d; n_2: %d" % (len(ivals), len(jvals)), fontsize=8, transform=ax.transAxes)

    plt.savefig("%s_ecdf.png" % propname)
    plt.close(fig)

    if (args.control):
        fig, ax = plt.subplots()

        mw_test_r = mannwhitneyu(rivals, rjvals)

        ri_sorted=np.sort( rivals )
        ri_yvals=np.arange(len(ri_sorted))/float(len(ri_sorted))
        plt.plot( ri_sorted, ri_yvals )

        rj_sorted=np.sort( rjvals )
        rj_yvals=np.arange(len(rj_sorted))/float(len(rj_sorted))
        plt.plot( rj_sorted, rj_yvals )

        plt.legend([propname + "_rand1", propname + "_rand2"], fontsize=8)
        plt.text(.7, 0.1, "M-W U p_ctrl = %.3g" % (mw_test_r[1]), fontsize=8, transform=ax.transAxes)

        propfile.write("set1,set2,rand1,rand2,ecdf_1_x,ecdf_1_y,ecdf_2_x,ecdf_2_y,ecdf_r1_x,ecdf_r1_y,ecdf_r2_x,ecdf_r2_y\n")
        propfile.write("\n".join([mkline(i) for i in izip_longest(ivals, jvals, rivals, rjvals, i_sorted, i_yvals, j_sorted, j_yvals, ri_sorted, ri_yvals, rj_sorted, rj_yvals)]))

        plt.savefig("%s_ctrl_ecdf.png" % propname)
        plt.close(fig)

    else:
        propfile.write("set1,set2,ecdf_1_x,ecdf_1_y,ecdf_2_x,ecdf_2_y\n")
        propfile.write("\n".join([mkline(i) for i in izip_longest(ivals, jvals, i_sorted, i_yvals, j_sorted, j_yvals)]))

    propfile.close()

# --- output tx-specific values for each set ---

set1_outfile = safe_open_file("_".join(args.set1[1::2]) + ".csv")
set2_outfile = safe_open_file("_".join(args.set2[1::2]) + ".csv")

header = ["txid"]

for propfile, propdict in txome_props.iteritems():
    for k in range(len(propdict["header"])):
        propname = propfile + "_" + propdict["header"][k]
        header.append(propname)

set1_outfile.write(",".join(header)  + "\n")
set1_outfile.write("\n".join([",".join([txid] + [tx_to_gene[txid]] + tx_list_set1[txid]) for txid in tx_list_set1.iterkeys()]))

set2_outfile.write(",".join(header)  + "\n")
set2_outfile.write("\n".join([",".join([txid] + [tx_to_gene[txid]] + tx_list_set2[txid]) for txid in tx_list_set2.iterkeys()]))

