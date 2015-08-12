#!/usr/bin/env python

# fpkm_to_tpm.py

# convert a file with an arbitrary number of columns to tpm using the following formula:

#   TPM_i = FPKM_i * 1e6 / sum(FPKM_g for all genes g)

# formula based on a Lynch paper http://lynchlab.uchicago.edu/publications/Wagner,%20Kin,%20and%20Lynch%20%282012%29.pdf



import sys, argparse, csv
from collections import defaultdict

parser = argparse.ArgumentParser(description="Convert FPKM/RPKM to TPM")
parser.add_argument("-i", "--input", help="File containing identifiers to use for the merge", required=True)
parser.add_argument("-t", "--separator", help="Field separator (default comma; \"tab\" for tabs; \"space\" for whitespace", default=",")
parser.add_argument("-o", "--output", nargs="?", help="File to output to (default stdout)", type=argparse.FileType('w'),
                    default=sys.stdout)
parser.add_argument("--ignore", help="Number of columns to ignore (one-based; 1 ignores the first column)", default=1, type=int)
parser.add_argument("--filter", help="Filter genes with TPM below arg", default=-1, type=float)
parser.add_argument("-u", "--unique", help="Only output lines with unique entries in column 1", action="store_true")

args = parser.parse_args()

if args.separator == "tab":
    delim="\t"
elif args.separator == "space":
    delim=" "
else:
    delim=args.separator

print "Converting FPKM values in %s to TPM" % args.input

infile = open(args.input, "r")
inreader = csv.reader(infile, quoting=csv.QUOTE_NONE, delimiter=delim)
out = csv.writer(args.output, delimiter=delim)

sums = defaultdict(float) 

header = inreader.next() 
print "Interpreting this line as a header:\n%s" % "\t".join(header)

for line in inreader: 
    for i in range(args.ignore,len(line)):
        sums[i] += float(line[i])

print sums

infile.seek(0)
inreader.next() 

out.writerow(header)

lastid = None

for line in inreader: 
    if (sum([float(i) for i in line[args.ignore:]])/float(len(line[args.ignore:]) > args.filter)):
        if (args.unique and line[0] != lastid):
            #print "line %s last %s" % (line[0], lastid)
            out.writerow(line[:args.ignore] + ["%6f" % (float(line[i]) * 1e6 / sums[i]) for i in range(args.ignore,len(line))])
            lastid = line[0]
        elif not args.unique:
            out.writerow(line[:args.ignore] + ["%6f" % (float(line[i]) * 1e6 / sums[i]) for i in range(args.ignore,len(line))])
