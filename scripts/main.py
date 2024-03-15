from Bio import SeqIO
from collections import Counter
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd

import kmers

#use RAINBOW as testing

reads = [s for s in SeqIO.parse("READS.fasta",'fasta')]
initial_queries = [r for r in SeqIO.parse("QUERY.fasta",'fasta')]

kmer_size = 3
# reads_dict = {}
# query_dict = {}

# make initial read k-mers
# for read in reads:
#     read_str = str(read.seq)
#     initial_r_kmers = kmers.kmers.make_kmers(read_str, kmer_size)
#     reads_dict[read.id] = initial_r_kmers
#
# for query in initial_queries:
#     query_str = str(query.seq)
#     initial_q_kmers = kmers.kmers.make_kmers(query_str, kmer_size)
#     query_dict[query.id] = initial_q_kmers

reads_dict = {1: ["NBO", "BOW"],
              2: ["UQH", "QHC"],
              3: ["AIN", "INB"],
              5: ["INB", "NBO"]
              }

query_dict = {6: ["RAI", "AIN"]
              }

# #match query kmers to read kmers
# #at the end add query to relavent reads dict
relavent_reads = {}

for qid, qkmers in query_dict.items():
    kmers.kmers.match_kmers(qkmers, reads_dict, relavent_reads)
    relavent_reads[qid] = qkmers

i = 0
while i < len(relavent_reads):
    rrkmers = list(relavent_reads.values())[i]
    kmers.kmers.match_kmers(rrkmers, reads_dict, relavent_reads)
    i += 1


print(relavent_reads)




### INITIAL DATA ANALYSIS ###
# seq_lengths = []
# for sequence in sequences:
#     seq_lengths.append(len(sequence.seq))

# average_length = sum(seq_lengths) / len()
# print(average_length)
# total_num_sequences = len(seq_lengths)



