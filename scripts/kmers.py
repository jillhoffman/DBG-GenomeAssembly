class kmers:
    # def __init__(self, read, kmer_size):
    #     self.read = read
    #     self.kmer_size = kmer_size

    def make_kmers(read, kmer_size):
        tot_length = len(read)
        all_kmers = []

        for i in range(tot_length):
            kmer = read[i:kmer_size]
            kmer_size += 1
            all_kmers.append(kmer)
            if kmer_size == tot_length + 1:
                break
        return all_kmers

    def match_kmers(kmers, reads_dict, relavent_reads):
        #relavent_reads[id] = kmers
        for kmer in kmers:
            for rid, rkmers in reads_dict.items():
                if kmer in rkmers:
                    if rid not in relavent_reads.keys():
                        relavent_reads[rid] = set(rkmers)