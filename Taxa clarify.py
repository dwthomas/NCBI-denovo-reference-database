#!/usr/bin/env python
### https://github.com/dwthomas/MiFish-Reference-Database
import pandas as pd
from Bio import SeqIO, pairwise2
from collections import Counter
import sys
import matplotlib.pyplot as plt
import progressbar
import os
import io
import subprocess
import gzip
import glob
import threading
from multiprocessing.pool import ThreadPool

### eg: python taxa_clarify_klymus_all.py \
#   klymus_full.tax \
#   klymus_full.fasta \
#   length_filtered_klymus_full \
#   merged_sequences/dna-sequences.fasta \
#   32 (number of threads to use)

taxa_in = sys.argv[1]
fasta_in = sys.argv[2]
tag = sys.argv[3]
merged_fasta = sys.argv[4]
threads = sys.argv[5]

fasta_out = tag + '.trimmed.fasta'
taxout = tag + ".tax"
final_fasta = tag + "trimmed_final.fasta"

header = "query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score".split(", ")

def run_blast(query):
    id = query.id
    folder_path = "blastout/" + id + "/"
    subject_file_path = folder_path + "subject.fasta"
    output_file_path = folder_path + "blastout"
    query_file_path = merged_fasta
    if os.path.exists(output_file_path+".gz"):
        try:
            with gzip.open(output_file_path+".gz", "rt") as btxt:
                best_hit = pd.read_csv(btxt, sep = "\t", names= header, index_col = 1)
        except:
            os.makedirs(folder_path, exist_ok = True)
            SeqIO.write([query], subject_file_path, "fasta")
            command = [
                "blastn",
                "-query", query_file_path,
                "-subject", subject_file_path,
                "-strand", "plus",
                "-task", "blastn",
                "-outfmt", "6",
                "-out", output_file_path,
                "-perc_identity", "0.6",
                "-qcov_hsp_perc", "0.6",
                "-max_target_seqs", "5",
            ]
            subprocess.run(command)
            best_hit = read_and_compress(output_file_path)
    else:
        os.makedirs(folder_path, exist_ok = True)
        SeqIO.write([query], subject_file_path, "fasta")
        command = [
            "blastn",
            "-query", query_file_path,
            "-subject", subject_file_path,
            "-strand", "plus",
            "-task", "blastn",
            "-outfmt", "6",
            "-out", output_file_path,
            "-perc_identity", "0.6",
            "-qcov_hsp_perc", "0.6",
            "-max_target_seqs", "5",
        ]
        subprocess.run(command)
        best_hit = read_and_compress(output_file_path)
    if len(best_hit) > 0:
        best_hit = best_hit.sort_values("evalue").iloc[0]
        return query[best_hit["s. start"]:best_hit["s. end"]], best_hit["bit score"]
    else:
        return None, 0

def compress(name, data):
    t = threading.Thread(target=os.remove, args=(name,))
    t.start()
    with gzip.open(name + ".gz", "wt") as bgz:
        bgz.write(data)
    t.join()

def read_and_compress(bout, r=True):
    with open(bout) as btxt:
        btxt = btxt.read()
    t = threading.Thread(target=compress, args=(bout, btxt))
    t.start()
    if r:
        pdout = pd.read_csv(io.StringIO(btxt), sep = "\t", names= header, index_col = 1)
    t.join()
    if r:
        return pdout

def r_and_c(bout):
    read_and_compress(bout, r=False)

taxa = pd.read_csv(taxa_in, sep = "\t", index_col = 0, names = ["Taxon"])

to_drop = []
to_keep = []
for seq in SeqIO.parse(fasta_in, "fasta"):
    l = len(seq)
    if l < 100:
        to_drop.append(seq.id)
    else:
        to_keep.append(seq)

SeqIO.write(to_keep, fasta_out, "fasta")

#causing error. Devin said to comment it out
taxa = taxa.drop(to_drop)
taxa.to_csv(taxout, sep = "\t", header=False)

outseqs = []
oute = []

seqs = list(SeqIO.parse(fasta_out, "fasta"))
with ThreadPool(threads) as pool:
    out = pool.map(run_blast, seqs)

for o in out:
    if hasattr(o[0], 'id'):
        outseqs.append(o[0])
        oute.append(o[1])

SeqIO.write(outseqs, final_fasta, "fasta")
