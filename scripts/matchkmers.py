class MatchKmers:

    def __init__(self, kmers, reads_dict):
        self.kmers = kmers
        self.reads_dict = reads_dict


    def match_kmers(self, relavent_reads, rr_overlap, q):

        for rid, rkmers in self.reads_dict.items():
            tmp = []

            for rkmer in rkmers:

                if rkmer in self.kmers:
                    tmp.append(rkmer)

            similarity = len(tmp)/len(rkmers)

            if 0.1 < similarity <= 0.9:
                relavent_reads[rid] = rkmers

                if q=="T":
                    rr_overlap[rid] = similarity

                else:
                    rr_overlap[rid] = 0

        return relavent_reads
