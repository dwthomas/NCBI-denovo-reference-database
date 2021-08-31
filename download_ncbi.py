#!/usr/bin/env python3
from Bio import Entrez  
import progressbar
import sys
import pickle


Entrez.email = sys.argv[1]
Entrez.api_key = sys.argv[2]
query = sys.argv[3]



def get_taxid(result):
    return result['LinkSetDb'][0]["Link"][0]["Id"]

def format_taxa(text):
    try:
        lineage = text["Lineage"]
    except:
        lineage = ""
    try:
        sciname = text["ScientificName"]
    except:
        sciname = ""
    try:
        cname = text['OtherNames']['GenbankCommonName']
    except:
        cname = ""
    return "; ".join((lineage, sciname, cname))

def get_accesions(query):
    print("Searching for matching accessions")
    retmax = 10000
    handle = Entrez.esearch(db="nucleotide",
                            retmax=retmax,
                            term=query,usehistory = "y",
                            idtype="acc")

    record = Entrez.read(handle)
    handle.close()

    total = int(record["Count"])
    accesions = set(record["IdList"])
    retstart = retmax
    webenv = record["WebEnv"]

    print("Found " + str(total) + " to download.", flush = True)

    with progressbar.ProgressBar(max_value=total) as bar:
        if retstart <= total:
            bar.update(retstart)
        else:
            bar.update(total)    
        while retstart <= total:
            handle = Entrez.esearch(db="nucleotide", 
                                    retstart = retstart,
                                    retmax=retmax, 
                                    term="(((12S OR MiFish) OR mitochondrion) NOT chromosome) NOT shotgun",
                                    usehistory = "y",
                                    webenv = webenv,
                                    idtype="acc")
            record = Entrez.read(handle)
            handle.close()
            accesions = accesions.union(set(record["IdList"]))
            retstart += retmax
            if retstart <= total:
                bar.update(retstart)
            else:
                bar.update(total)
    return list(accesions), webenv

def get_taxids(accesions, webenv=None):
    print("Fetching corresponding taxonomy ids", flush = True)
    taxids = []
    retmax = 1000
    total = len(accesions)
    with progressbar.ProgressBar(max_value=total) as bar:
        retstart = 0
        if retstart <= total:
            bar.update(retstart)
        else:
            bar.update(total)
        while retstart <= len(accesions):
            handle = Entrez.elink(db = "taxonomy",
                                  dbfrom = "nucleotide",
                                  id = accesions[retstart:min(len(accesions),retmax + retstart)],
                                  webenv = webenv)
            results = Entrez.read(handle)
            handle.close()
            retstart += retmax
            taxids += list(map(get_taxid, results))
            if retstart <= total:
                bar.update(retstart)
            else:
                bar.update(total)
    return taxids

def fetch_taxa(tids, webenv=None):
    print("Downloading taxonomic information.", flush = True)
    taxa_info = {}
    retmax = 1000
    total = len(tids)
    with progressbar.ProgressBar(max_value=total) as bar:
        retstart = 0
        if retstart <= total:
            bar.update(retstart)
        else:
            bar.update(total)
        while retstart <= total:
            handle = Entrez.efetch(db = "taxonomy", webenv = webenv, id = tids[retstart:min(total,retmax + retstart)])
            result = Entrez.read(handle)
            handle.close()
            for i in result:
                taxa_info[i["TaxId"]] = format_taxa(i)
            retstart += retmax
            if retstart <= total:
                bar.update(retstart)
            else:
                bar.update(total)
    return taxa_info

def fetch_seq(accesion, webenv):
    print("Downloading sequences", flush = True)
    seqs = []
    taxids = []
    retmax = 100
    total = len(accesion)
    with progressbar.ProgressBar(max_value=total) as bar:
        retstart = 0
        if retstart <= total:
            bar.update(retstart)
        else:
            bar.update(total)
        while retstart <= total:
            accbatch = accesion[retstart:min(total,retmax + retstart)]
            handle = Entrez.efetch(db = "nucleotide", webenv = webenv, rettype="fasta", retmode = "xml", id = accbatch)
            results = Entrez.read(handle)
            handle.close()
            seqs += list(map(lambda x: x["TSeq_sequence"], results))
            taxids += list(map(lambda x: x['TSeq_taxid'], results))

            with open("12Snmito.fasta" , "at") as seqout:
                with open("12Snmito.tax", "at") as taxout:
                    for i in range(len(accbatch)):
                        print(">" + accbatch[i], file = seqout)
                        print(seqs[retstart+i], file = seqout)
            
                        print(accbatch[i], end = "\t", file = taxout)
                        print(taxids[retstart+i], file = taxout)

            retstart += retmax
            if retstart <= total:
                bar.update(retstart)
            else:
                bar.update(total)
    return taxids, seqs

#accesions, webenv = get_accesions(query)
already_got = []
taxids = []
with open("12Snmito.tax") as taxfile:
    for line in taxfile:
        already_got.append(line.split()[0])
        taxids.append(line.split()[1])
#accesions = list(set(accesions) - set(already_got))
#taxids, seqs = fetch_seq(accesions, webenv)
tids = list(set(taxids))
taxa_info = {}
print(len(taxa_info), len(tids))
taxa_info.update(fetch_taxa(tids))
tids = list(set(tids) - set(taxa_info.keys()))
#print(tids)
with open('taxa_info.pickle', 'wb') as handle:
        pickle.dump(taxa_info, handle)
with open("12SnMito_full.tax", "wt") as ft:
    for i in range(len(taxids)):
        if taxids[i] in taxa_info: 
            print(already_got[i] + "\t" + taxa_info[taxids[i]], file = ft)
        else:
            try:
                print(already_got[i], taxids[i], i in tids)
                second_tid = get_taxids([already_got[i]])
                print(second_tid)
                gst = fetch_taxa(second_tid)
                print(gst)
                print(already_got[i] + "\t" + list(gst.values())[0], file = ft)
            except:
                print(already_got[i] + "\t" + already_got[i], file = ft)
                


print("Creating files")
