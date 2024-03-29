from Bio import SeqIO
import makekmers
import sys
import matchkmers
import order
import time

def main():
    make_kmer_st = time.time()

    # reads = [s for s in SeqIO.parse("./../example_data/READS.01.fasta", 'fasta')]
    # initial_queries = [r for r in SeqIO.parse("./../example_data/QUERY.fasta", 'fasta')]
    # kmer_size = 4

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
        rkmercuts = makekmers.MakeKmers(read_str, kmer_size)
        initial_r_kmers = rkmercuts.make_kmers()
        reads_dict[read.id] = initial_r_kmers


    for query in initial_queries:
        query_str = str(query.seq)
        qkmercuts = makekmers.MakeKmers(query_str, kmer_size)
        initial_q_kmers = qkmercuts.make_kmers()
        query_dict[query.id] = initial_q_kmers

    num_reads = len(reads_dict)
    print(f"matching {num_reads} reads to initial query...")
    relavent_reads = {}
    rr_overlap = {}

    for qid, qkmers in query_dict.items():

        matching = matchkmers.MatchKmers(qkmers, reads_dict)
        matching.match_query_kmers(relavent_reads, rr_overlap)
        relavent_reads[qid] = qkmers


    query_rr = len(relavent_reads) - 1

    if len(rr_overlap.values()) != 0:
        max_overlap = max(rr_overlap.values())
        id_max_overlap = list(filter(lambda x: rr_overlap[x] == max_overlap, rr_overlap))

        print(f"matched {query_rr} reads, the highest at ({max_overlap*100}%) kmer overlap")
        print(f"checking for more relevent reads...")

        for id in id_max_overlap:
            rrkmers = relavent_reads[id]
            rrmatching = matchkmers.MatchKmers(rrkmers, reads_dict)
            rrmatching.match_lr_kmers(relavent_reads)

        print(f"found {len(relavent_reads) - query_rr - 1} additional reads overlapped")
    else:
        print(f"no reads were matched, try a smaller k-mer size")

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
    # for uniqueKmer in unique_kmers:
    #     km1cuts = makekmers.MakeKmers(uniqueKmer, km1_size)
    #     km1 = km1cuts.make_kmers()
    #     km1_dict[uniqueKmer] = km1
    #
    # total_kmers = len(km1_dict)
    # ordered_kmers = {}
    # current_order = 1
    #
    #
    # others = order.Order(km1_dict)
    #
    # for kmer, km1ers in km1_dict.items():
    #     others.otherKm1ers(kmer, km1ers, ordered_kmers, current_order)
    #
    #
    # while current_order <= total_kmers:
    #     c_kmer = ordered_kmers[current_order]
    #     c_km1ers = km1_dict[c_kmer]
    #     next_order = current_order + 1
    #     others.orderKmers(c_km1ers, c_kmer, ordered_kmers, next_order)
    #     current_order += 1

if __name__ == '__main__':
    main()




