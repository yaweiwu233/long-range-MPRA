import sys
import os
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
import re
from optparse import OptionParser
from collections import defaultdict
import pandas as pd
import numpy as np

parser = OptionParser()
parser.add_option("-f", "--file",dest="f",help="fastq file name")
parser.add_option("-e", "--eBC", dest="e", help="enhancer barcode file directory")
parser.add_option("-d", "--dBC", dest="d", help="distance barcode file directory")
parser.add_option("-o", "--out", dest="o",help="output file directory", default="./")
parser.add_option("--eBC-pattern", dest="pe", help="search pattern for enhancer barcode")
parser.add_option("--dBC-pattern", dest="pd", help="search pattern for distance barcode")

options, args = parser.parse_args()

filename = options.f
basename = os.path.basename(filename)
eBC_dir = options.e
dBC_dir = options.d
outdir = options.o
eBC_pattern = re.compile(options.pe)
dBC_pattern = re.compile(options.pd)

print("sample_ID:", basename)
print("eBC_pattern:", eBC_pattern)
print("dBC_pattern:", dBC_pattern)

def read_fastq(fastq_file):

    current_record = {}

    for name, seq, crap, quality in zip(*[iter(fastq_file)]*4):
        current_record['name'] = name.strip('\n')
        current_record['seq'] = seq.strip('\n')
        current_record['quality'] = quality.strip('\n')

        yield current_record

def get_fastq_id(fastq_name):

    return fastq_name.split(' ')[0]

def rc(seq):
    complement = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    return ''.join([complement[x] for x in seq[::-1]])
    

def fuzzy_bc(query, barcodes, threshold):
    for x in barcodes:
        if x in query:
            return x
        score = fuzz.partial_ratio(x, query)
        if score >= threshold:
            return x
    return False

# Read eBCs and dBCs into dictionaries
f2 = open(dBC_dir, "r")

eBCs = {}
dBCs = {}

eBC_file = pd.read_csv(eBC_dir, sep="\t")
eBC_names = eBC_file["name"].to_list()
eBC_seqs = eBC_file["eBC"].to_list()
eBCs = {x:y for x,y in zip(eBC_seqs, eBC_names)}
for line in f2:
    x = line.strip().split("\t")
    dBCs[x[1]] = x[0]
    
f2.close()

# Start process fastq file
R1_file = gzip.open(filename+"_R1.fastq.gz", "rt")
R2_file = gzip.open(filename+"_R2.fastq.gz", "rt")

R1_reader = read_fastq(R1_file)
R2_reader = read_fastq(R2_file)

trios = defaultdict(int)
names = []
locations = []
integrations = []

i = 0
j = 0
k = 0

for R1_record, R2_record in zip(R1_reader, R2_reader):
    i += 1
    R1_seq = R1_record["seq"]
    R2_seq = R2_record["seq"]
    # parse cBC and rBC
    # Check read length, discarded those are shorter than 113
    if len(R1_seq) < 90 and len(R2_seq) < 50:
        continue
    # search for CRE barcodes and genome barcodes
    m = eBC_pattern.search(R1_seq)
    n = dBC_pattern.search(R2_seq)
    if (m == None) or (n == None):
        continue
    j += 1
    eBC = m.group(1)
    rBC = m.group(2)
    dBC = n.group(1)
    if (eBC in eBCs) and (dBC in dBCs):
        trios[(eBCs[eBC], dBCs[dBC], rBC)] += 1
        k += 1
print("total reads:", i)
print("reads with pattern", j)
print("reads with valid eBC and dBC:", k, k/i)

summary = pd.DataFrame(list(trios.keys()), columns = ["CRE", "LP", "rBC"])
summary["count"] = list(trios.values())
summary.to_csv(outdir+"/"+basename+"_barcode_summary.txt", sep = "\t", index = False)
integrations = summary.groupby(["CRE", "LP"])["rBC"].nunique()
print("Number of integrations:")
print(integrations)
print("Number of reads:")
reads = summary.groupby(["CRE", "LP"])["count"].sum()
print(reads)

# Plot read distribution across different cBC-gBC pair
sns.histplot(x=summary.groupby(["CRE", "LP"]).sum()["count"].to_list(), log_scale=True)
plt.savefig(outdir+"/"+basename+"_eBC-dBC_count_distribution.png", dpi=300)

