class MatchKmers:

    def __init__(self, kmers, reads_dict):
        self.kmers = kmers
        self.reads_dict = reads_dict


    def match_query_kmers(self, relavent_reads, rr_overlap):
        for rid, rkmers in self.reads_dict.items():
            tmp = []
            for qkmer in self.kmers:
                if qkmer in rkmers:
                    tmp.append(qkmer)

            similarity = len(tmp)/len(self.kmers)
            if similarity >= 0.5:
                relavent_reads[rid] = rkmers
                rr_overlap[rid] = similarity

        return relavent_reads


    def match_lr_kmers(self, relevant_reads):

        for rid, rkmers in self.reads_dict.items():
            tmp = []
            for qkmer in self.kmers:
                if qkmer in rkmers:
                    tmp.append(qkmer)

            similarity = len(tmp) / len(self.kmers)
            if similarity >= 0.5:
                relevant_reads[rid] = set(rkmers)

        return relevant_reads