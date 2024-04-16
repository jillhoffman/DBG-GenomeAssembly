class MakeKmers:

    def __init__(self, read, kmer_size):
        self.read = read
        self.kmer_size = kmer_size

    def make_kmers(self):

        tot_length = len(self.read)
        all_kmers = []

        for i in range(tot_length):
            kmer = self.read[i:self.kmer_size]
            self.kmer_size += 1
            all_kmers.append(kmer)

            if self.kmer_size == tot_length + 1:
                break

        return all_kmers


    def make_km1ers(self, qkmers):

            tmp = ()
            tot_length = len(self.read)

            for i in range(tot_length):
                Kmer=self.read

                if Kmer not in qkmers and Kmer != qkmers[-1:]:

                    km1er = Kmer[i:self.kmer_size]
                    self.kmer_size += 1
                    tmp = tmp + (km1er,)

                if self.kmer_size == tot_length + 1:
                    break

            return tmp










