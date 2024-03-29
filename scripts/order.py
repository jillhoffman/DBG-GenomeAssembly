class Order():

    def __init__(self, km1_dict):
        self.km1_dict = km1_dict


    def otherKm1ers(self, kmer, km1ers, ordered_kmers, cOrder):

        firstKm1ers = []
        secondKm1ers = []

        for kmer2, km1ers2 in self.km1_dict.items():
            if kmer != kmer2:
                firstKm1ers.append(km1ers2[0])
                secondKm1ers.append(km1ers2[1])

        if km1ers[1] in firstKm1ers and km1ers[0] not in secondKm1ers:
            ordered_kmers[cOrder] = kmer

        return ordered_kmers

    def orderKmers(self, current_km1ers, current_kmer, ordered_kmers, nOrder):

        for kmer, k1mers in self.km1_dict.items():
            if current_km1ers[1] == k1mers[0] and kmer != current_kmer:
                if nOrder not in list(ordered_kmers.keys()):
                    ordered_kmers[nOrder] = kmer
                else:
                    ordered_kmers[nOrder] = [ordered_kmers[nOrder], kmer]