#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO, pairwise2
from collections import Counter
import matplotlib.pyplot as plt
import progressbar
import os
import io
import subprocess
import gzip
import glob
import threading
from multiprocessing.pool import ThreadPool

taxa = pd.read_csv("./12SnMito_full.tax", sep = "\t", index_col = 0, names = ["Taxon"])

to_drop = []
to_keep = []
for seq in SeqIO.parse("12Snmito.fasta", "fasta"):
    l = len(seq)
    if l < 100:
        to_drop.append(seq.id)
    else:
        to_keep.append(seq)
SeqIO.write(to_keep, "length_filtered_12SnMito.fasta", "fasta")

taxa = taxa.drop(to_drop)
taxa.to_csv("length_filtered_12SnMito.tax", sep = "\t")


header = "query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score".split(", ")
def run_blast(query):
    id = query.id
    folder_path = "blastout/" + id + "/"
    subject_file_path = folder_path + "subject.fasta"
    output_file_path = folder_path + "blastout"
    query_file_path = "./reference-seqs.fasta"
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

outseqs = []
oute = []
seqs = list(SeqIO.parse("length_filtered_12SnMito.fasta", "fasta"))
with ThreadPool(72) as pool:
    out = pool.map(run_blast, seqs)
for o in out:
    outseqs.append(o[0])
    oute.append(o[1])
