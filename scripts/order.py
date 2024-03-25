class Order():
    def all_other_km1ers(self, km1_dict, kmer):

        all_others = []
        for kmer2, km1ers2 in km1_dict.items():
            if kmer != kmer2:
                for k1mer2 in km1ers2:
                    all_others.append(k1mer2)

        return all_others

    def order_kmers(self, km1_dict, current_km1ers, current_kmer, ordered_kmers, order):

        for kmer, k1mers in km1_dict.items():
            if current_km1ers[1] in k1mers and kmer != current_kmer:
                ordered_kmers[order] = kmer