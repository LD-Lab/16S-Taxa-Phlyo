#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Thanks for the original author of this scripts:
# @author: S. Kim


import os
import re
import argparse
import numpy as np
# Biopython
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


""" Initialize environment """

#WD = os.path.expanduser("~/Projects/David Lab/Ind study projects/Seq trimming/Data/")
# Order by fwd/rev pairs. Even indices are fwd, odd indices are rev.
# Note that all rev sequences will be reverse-complemented to allow merging,
# which means all indices for rev reads will be mirrored.

IMPORTED_FP = "seqs.fasta"
TRIMMED_FP = "trimmed_seqs.fasta"
MERGED_FP = "merged_seqs.fasta"
UNMERGED_FP = "unmerged_seqs.fasta"
BLAST_FP = "blast_results.txt"


# Import all paired forward/reverse sequences.
# Each dict key is paired with a 2-tuple, corresponding to fwd/rev.
# wd: Working directory, where reads are stored
# names: Filenames
# return: [{names, reads, seqs, quals, trimmed_seqs, trim_indices, merged_seq}]
def import_seqs(filenames, Flist, Rlist):
    sequences = []
    for i in range(len(filenames)):
        names = (filenames[i]+"_F", filenames[i]+"_R")
        reads = (SeqIO.read(Flist[i], "abi"),
                 SeqIO.read(Rlist[i], "abi"))
        seqs = (str(reads[0].seq), rev_cmp(str(reads[1].seq)))
        quals = ([x for x in reads[0].letter_annotations["phred_quality"]],
                 [x for x in reads[1].letter_annotations["phred_quality"]][::-1])
        sequences.append({"names": names, "seqs": seqs, "quals": quals})
    return sequences


# Find the reverse complement of a DNA sequence.
# seq: Sequence
# return: Rev cmp sequence
DNA_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def rev_cmp(seq):
    revcmp_seq = ""
    for c in reversed(seq):
        revcmp_seq += DNA_COMP[c]
    return revcmp_seq


# Trim sequence by quality scores, using Kadane's algorithm.
# seq: Raw sequence
# quality_scores: Phred quality scores for seq
# (Optional) name: Name of seq
# (Optional) fwd: Forward/reverse, used to position trimmed seq in original seq
# (Optional) cutoff: Lowest acceptable quality score, used to normalize scores around 0
# (Optional) low_penalty: Penalty per nt of including a low-quality base
# return: Trimmed seq, position in original seq
def trim_seq(seq, quality_scores, name=None, fwd=True, cutoff=None, low_penalty=None):
    if name is not None:
        print(name)
    if cutoff is None:
        # mean log, better cutoff metric needed
        cutoff = np.average(quality_scores) * 0.8
    if low_penalty is None:
        low_penalty = cutoff * -3  # arbitrary gap penalty, seems to work
    scores = [(x-cutoff if x >= cutoff else x-cutoff+low_penalty)
              for x in quality_scores]

    # Kadane's algorithm
    max_to_here = 0
    max_so_far = float("-inf")
    ind_to_here = [0, -1]
    ind_so_far = [0, -1]
    for i in range(len(seq)):
        max_to_here += scores[i]
        ind_to_here[1] += 1
        if max_to_here > max_so_far:
            max_so_far = max_to_here
            ind_so_far = list(ind_to_here)
        if max_to_here < 0:
            max_to_here = 0
            ind_to_here = [i+1, i]

    if max_to_here <= 0 and max_so_far <= 0:
        print("Entire sequence was low quality. Try a lower cutoff.\n")
        return ("", "")
    max_ind = ind_to_here if max_to_here > max_so_far else ind_so_far
    trimmed_seq = seq[max_ind[0]:max_ind[1]+1]
    trimmed_scores = quality_scores[max_ind[0]:max_ind[1]+1]
    print("Trimmed seq: {} nt ({}-{}, {:.0f}%)".format(
        len(trimmed_seq),
        max_ind[0]+1 if fwd else len(seq)-max_ind[1],
        max_ind[1]+1 if fwd else len(seq)-max_ind[0],
        100*len(trimmed_seq)/len(seq)))
    print("Avg score: {:.2f} (std {:.2f}, rng {}-{}, >={:.0f}, {:.0f})\n".format(
        np.average(trimmed_scores), np.std(trimmed_scores),
        min(trimmed_scores), max(trimmed_scores), cutoff, low_penalty))
    return trimmed_seq, max_ind


# Merge forward and reverse reads, using the Smith-Waterman algorithm.
# Different scores are assigned for matches, mismatches, and gaps in the alignment.
# fwd/rev: Forward/reverse trimmed reads
# (Optional) name_fwd/rev: Names of fwd/rev reads
# (Optional) ind_fwd/rev: Positions of trimmed reads in original seqs
# (Optional) seq_fwd/rev: Original fwd/rev sequences
# (Optional) match/mismatch/gap_penalty: Scores for alignment
# return: merged seq
def merge_seqs(fwd, rev, name_fwd=None, name_rev=None, ind_fwd=None, ind_rev=None,
               len_seq_rev=None, match=5, mismatch=-5, gap_penalty=-5):
    if name_fwd is not None and name_rev is not None:
        print(name_fwd + " + " + name_rev)

    # Build dictionary of scores for aligning bases
    dna_bases = ["A", "T", "C", "G", "N"]
    score_dict = {}
    for i in range(len(dna_bases)):
        for j in range(len(dna_bases)):
            score_dict[dna_bases[i]+dna_bases[j]
                       ] = match if i == j else mismatch
    score_dict["NN"] = mismatch  # no use matching unknown with unknown

    # Initialize dynamic programming and traceback tables
    LEFT, DIAG, UP = range(3)   # define pointers: 0==Left, 1==Diagonal, 2==Up
    dp_table = [[0]*(len(rev)+1) for _ in range(len(fwd)+1)]    # base case = 0
    I = len(dp_table)
    J = len(dp_table[0])
    tb_table = [[[] for _ in range(len(rev)+1)]
                for _ in range(len(fwd)+1)]  # no base case

    # Build tables
    max_align_score = 0
    max_back_index = (0, 0)
    for i in range(1, I):
        for j in range(1, J):
            scores = (gap_penalty + dp_table[i][j-1],  # left
                      score_dict[fwd[i-1]+rev[j-1]] + \
                      dp_table[i-1][j-1],   # diagonal
                      gap_penalty + dp_table[i-1][j],  # up
                      0)    # quit
            dp_table[i][j] = max_score = max(scores)
            if max_score != 0:
                for k in range(len(scores)):
                    if scores[k] == max_score:
                        tb_table[i][j].append(k)
                if max_score > max_align_score:
                    max_align_score = max_score
                    max_back_index = (i, j)

    # Traceback
    max_i, max_j = max_back_index
    max_front_index = (0, 0)
    fwd_align = rev_align = ""
    while len(tb_table[max_i][max_j]) > 0:
        max_front_index = (max_i, max_j)
        ptr = max(tb_table[max_i][max_j])
        if ptr != UP:
            if ptr == LEFT:
                fwd_align = "-" + fwd_align
            rev_align = rev[max_j-1] + rev_align
            max_j -= 1
        if ptr != LEFT:
            if ptr == UP:
                rev_align = "-" + rev_align
            fwd_align = fwd[max_i-1] + fwd_align
            max_i -= 1
    matches, indels = 0, 0
    for i in range(len(fwd_align)):
        if fwd_align[i] == rev_align[i]:
            matches += 1
        elif fwd_align[i] == "-" or rev_align[i] == "-":
            indels += 1

    # Merge forward, aligned overlap, and reverse, and compile stats
    merged_seq = ""
    merged_nts = ""
    discarded_nts = [0, 0]
    print_merge_stats = None not in (ind_fwd, ind_rev, len_seq_rev)
    if max_front_index[0] > max_front_index[1]:  # before overlap
        merged_seq += fwd[:max_front_index[0]]
        if print_merge_stats:
            merged_nts += "fwd {}-{}, ".format(
                ind_fwd[0]+1,
                ind_fwd[0]+max_front_index[0]-1)
        discarded_nts[1] += max_front_index[1]-1
    else:
        merged_seq += rev[:max_front_index[1]]
        if print_merge_stats:
            merged_nts += "rev {}-{}, ".format(
                len_seq_rev-ind_rev[0]+1,
                len_seq_rev-ind_rev[0]-max_front_index[1])
        discarded_nts[0] += max_front_index[0]-1
    merged_seq += fwd[max_front_index[0]:max_back_index[0]]  # overlap
    if print_merge_stats:
        merged_nts += "fwd {}-{}/rev {}-{}, ".format(
            ind_fwd[0]+max_front_index[0],
            ind_fwd[0]+max_back_index[0],
            len_seq_rev-ind_rev[0]-max_front_index[1]+1,
            len_seq_rev-ind_rev[0]-max_back_index[1]+1)
    if len(fwd)-max_back_index[0] < len(rev)-max_back_index[1]:  # after overlap
        merged_seq += rev[max_back_index[1]:]
        if print_merge_stats:
            merged_nts += "rev {}-{}".format(
                len_seq_rev-ind_rev[0]-max_back_index[1],
                len_seq_rev-ind_rev[1])
        discarded_nts[0] += len(fwd)-max_back_index[0]
    else:
        merged_seq += fwd[max_back_index[0]:]
        if print_merge_stats:
            merged_nts += "fwd {}-{}".format(
                ind_fwd[0]+max_back_index[0]+1,
                ind_fwd[1]+1)
        discarded_nts[1] += len(rev)-max_back_index[1]

    overlap_len = max_back_index[0]-max_front_index[0]+1
    if overlap_len <= 2*discarded_nts[0] or overlap_len <= 2*discarded_nts[1]:
        merged_seq = None
        print("Reads may not overlap or be too low quality. ({} nt, lost {} fwd/{} rev)\n".format(
            max_back_index[0]-max_front_index[0]+1, discarded_nts[0], discarded_nts[1]))
    else:
        print("Merged seq: {} nt ({})".format(len(merged_seq), merged_nts))
        print("Overlap: {} nt ({:.1f}% match, {} indels, lost {} fwd/{} rev)\n".format(
            max_back_index[0]-max_front_index[0]+1, 100*matches/len(fwd_align),
            indels, discarded_nts[0], discarded_nts[1]))
    return merged_seq


# BLAST batch of reads.
# fp: FASTA file path
# return: BLAST records
def blast_seqs(fp):
    # Submit BLAST query
    return list(NCBIXML.parse(NCBIWWW.qblast("blastn", "nt", open(fp).read())))


# Write FASTA file.
# sequences: Dictionary of reads
# filename: Filename
# field: original, trimmed, or merged
def write_fasta(sequences, filename, field=""):
    records = []
    for sq in sequences:
        if field == "original":
            for k in range(2):
                records.append(SeqRecord(
                    Seq(sq["seqs"][k], IUPAC.IUPACAmbiguousDNA()),
                    id=sq["names"][k], description=""))
        elif field == "trimmed" or sq["merged_seq"] is None:
            for k in range(2):
                records.append(SeqRecord(
                    Seq(sq["trimmed_seqs"][k], IUPAC.IUPACAmbiguousDNA()),
                    id=sq["names"][k], description="({})".format(
                        "unmerged" if field == "merged" else field)))
        elif field == "merged":
            records.append(SeqRecord(
                Seq(sq["merged_seq"], IUPAC.IUPACAmbiguousDNA()),
                id=", ".join(sq["names"]), description="(merged)"))
    SeqIO.write(records, filename, "fasta")


parser = argparse.ArgumentParser(description='Run Diamond')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
args = parser.parse_args()

inputDir = str(args.fileDir)
outputDir = os.path.abspath(args.OpDir)


""" Main """

# Import fwd/rev reads
#data = "/mnt/d/Lab/TaxaIdentification/data/16S-seq-result"
data = inputDir

FFileList = []
RFileList = []

# Find all forward sequence files and reverse sequence files
for file in os.listdir(data):
    filestr = re.sub("[-()[\]]", "_", file)
    filestr = re.sub("\.", "__", filestr)
    if re.search("__.*F__", filestr) and file.endswith(".ab1"):
        FFileList.append(file)
    elif re.search("__.*R__", filestr) and file.endswith(".ab1"):
        RFileList.append(file)

# Match the forward sequence files with reverse sequence files
filenames = []
Flist = []
Rlist = []
for FFile in FFileList:
    for RFile in RFileList:
        FFileStr = re.sub("[-()[\]]", "_", FFile)
        FFileStr = re.sub("\.", "__", FFileStr)
        RFileStr = re.sub("[-()[\]]", "_", RFile)
        RFileStr = re.sub("\.", "__", RFileStr)
        flag = 0
        filename = []
        filenamestr =[]
        for i in range(len(re.split("__+",FFileStr))-1):
            if re.split("__+",FFileStr)[i] == re.split("__+",RFileStr)[i] and re.split("__+",FFileStr)[i] != []:
                flag = 1
        if flag == 1:
            for i in range(len(re.split("_+",FFileStr))-1):
                if re.split("_+",FFileStr)[i] == re.split("_+",RFileStr)[i]:
                    filename.append(re.split("_+",FFileStr)[i])
            filenamestr = '_'.join(filename)
            filenames.append(filenamestr)
            Flist.append(os.path.join(data, FFile))
            Rlist.append(os.path.join(data, RFile))


# Check the "__***F__" pattern to get the strain name form the forward sequence file. If the strain name is too long, find the real strain name with a "DA***" pattern.

#for file in os.listdir(data):
#    if re.search("_[\[_].*F[\]_][._]", file) and file.endswith(".ab1"):
#        filename = file.split(re.findall(r"_[\[_].*F[\]_][._]", file)[0])[0]
#        if len(filename)>10 and "DA" in filename:
#            filename = re.findall("DA.\d+",filename)[0]
#        filenames.append(filename)
        #Fpath = os.path.join(data, file)
        #Rpath = os.path.join(data, file[:-10]+"_[1492R].ab1")

# Find and list the forward and reverse sequence file for each strain name.

#for i in range(len(filenames)):
#    for file in os.listdir(data):
#        if filenames[i] in file and file.endswith(".ab1"):
#            if re.search("_[\[_].*F[\]_][._]", file):
#                Flist.append(os.path.join(data, file))
#            elif re.search("_[\[_].*R[\]_][._]", file):
#                Rlist.append(os.path.join(data, file))

sequences = import_seqs(filenames, Flist, Rlist)
write_fasta(sequences, os.path.join(outputDir, IMPORTED_FP), "original")

# Trim fwd/rev reads
print("================== TRIMMING SEQUENCES ==================")
for _sq in sequences:
    _fwd_trim, _fwd_ind = trim_seq(_sq["seqs"][0], _sq["quals"][0], _sq["names"][0], True,
                                   30, -90)    # qual=30 == 99.9% accuracy
    _rev_trim, _rev_ind = trim_seq(_sq["seqs"][1], _sq["quals"][1], _sq["names"][1], False,
                                   30, -90)
    _sq["trimmed_seqs"] = (_fwd_trim, _rev_trim)
    _sq["trim_indices"] = (_fwd_ind,  _rev_ind)
write_fasta(sequences, os.path.join(outputDir, TRIMMED_FP), "trimmed")

# Merge reads
print("================== MERGING SEQUENCES ===================")
for _sq in sequences:
    _sq["merged_seq"] = merge_seqs(_sq["trimmed_seqs"][0], _sq["trimmed_seqs"][1],
                                   _sq["names"][0],        _sq["names"][1],
                                   _sq["trim_indices"][0], _sq["trim_indices"][1],
                                   len(_sq["seqs"][1]),
                                   2, -7, -7)   # arbitrary scores, seems to work
write_fasta(sequences, os.path.join(outputDir, MERGED_FP), "merged")

'''
# BLAST merged reads and parse output
print("================== BLASTING SEQUENCES ==================")
blast_records = blast_seqs(os.path.join(outputDir, MERGED_FP))
_k = 0.0
for _record in blast_records:
    if _record.query.endswith("(unmerged)"):
        if "blast" not in sequences[int(_k)]:
            sequences[int(_k)]["blast"] = []
            sequences[int(_k)]["blast_e"] = []
        if len(_record.alignments) > 0:
            sequences[int(_k)]["blast"].append(_record.descriptions[0].title)
            sequences[int(_k)]["blast_e"].append(_record.descriptions[0].e)
        else:
            sequences[int(_k)]["blast"].append("No matches")
            sequences[int(_k)]["blast_e"].append(float("-inf"))
        _k += 0.5
    else:
        sequences[int(_k)]["blast"] = _record.descriptions[0].title
        sequences[int(_k)]["blast_e"] = _record.descriptions[0].e
        _k += 1.0
for _sq in sequences:
    print("{} ({}merged)".format(
        ", ".join(_sq["names"]), "un" if _sq["merged_seq"] is None else ""))
    if type(_sq["blast"]) is list:
        print("BLAST fwd: {} (e={:.1f})".format(
            _sq["blast"][0], _sq["blast_e"][0]))
        print("BLAST rev: {} (e={:.1f})\n".format(
            _sq["blast"][1], _sq["blast_e"][1]))
    else:
        print("BLAST: {} (e={:.1f})\n".format(_sq["blast"], _sq["blast_e"]))
'''