from Bio import SeqIO
import makekmers_e
import sys
import matchkmers_e
import time

def main():
    make_kmer_st = time.time()

    fasta_queries = sys.argv[1]
    initial_queries = [r for r in SeqIO.parse(fasta_queries,'fasta')]
    fasta_reads = sys.argv[2]
    reads = [s for s in SeqIO.parse(fasta_reads,'fasta')]
    kmer_size = int(sys.argv[3])


    reads_dict = {}
    query_dict = {}
    print(f"creating {kmer_size}-mers...")

    for read in reads:
        read_str = str(read.seq)
        rkmercuts = makekmers_e.MakeKmers(read_str, kmer_size)
        initial_r_kmers = rkmercuts.make_kmers()
        reads_dict[read.id] = initial_r_kmers


    for query in initial_queries:
        query_str = str(query.seq)
        qkmercuts = makekmers_e.MakeKmers(query_str, kmer_size)
        initial_q_kmers = qkmercuts.make_kmers()
        query_dict[query.id] = initial_q_kmers

    num_reads = len(reads_dict)
    print(f"matching {num_reads} reads to initial query...")
    relavent_reads = {}

    for qid, qkmers in query_dict.items():

        matching = matchkmers_e.MatchKmers(qkmers, reads_dict)
        matching.match_query_kmers(relavent_reads)
        relavent_reads[qid] = qkmers

    query_rr = len(relavent_reads) - 1

    i = 0
    while i < len(relavent_reads):
        rrkmers = list(relavent_reads.values())[i]
        rrmatching = matchkmers_e.MatchKmers(rrkmers, reads_dict)
        rrmatching.match_lr_kmers(relavent_reads)
        i += 1

    #print(f"found {len(relavent_reads) - query_rr - 1} additional reads overlapped")


    print(relavent_reads)
    print("time elapsed: %s seconds" % (round(time.time() - make_kmer_st)))

    # km1_size = kmer_size - 1
    # unique_kmers = []
    # km1_dict = {}
    #
    # for rr_kmers in relavent_reads.values():
    #     for rr_kmer in rr_kmers:
    #         if rr_kmer not in unique_kmers:
    #             unique_kmers.append(rr_kmer)
    #
    # for unique_kmer in unique_kmers:
    #     km1 = kmers.kmers.make_kmers(unique_kmer, km1_size)
    #     km1_dict[unique_kmer] = km1

    # total_kmers = len(km1_dict)
    # ordered_kmers = {}
    # current_order = 1
    #
    # for kmer, km1ers in km1_dict.items():
    #     all_others = kmers.kmers.all_other_km1ers(km1_dict, kmer)
    #     if km1ers[1] in all_others and km1ers[0] not in all_others:
    #         ordered_kmers[current_order] = kmer
    #
    # while current_order <= total_kmers:
    #     c_kmer = ordered_kmers[current_order]
    #     c_km1ers = km1_dict[c_kmer]
    #     next_order = current_order + 1
    #     kmers.kmers.order_kmers(km1_dict, c_km1ers, c_kmer, ordered_kmers, next_order)
    #     current_order += 1

if __name__ == '__main__':
    main()




