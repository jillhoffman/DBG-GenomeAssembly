class MatchKmers:

    def __init__(self, kmers, reads_dict):
        self.kmers = kmers
        self.reads_dict = reads_dict


    def match_query_kmers(self, relavent_reads):
        for rid, rkmers in self.reads_dict.items():
            for qkmer in self.kmers:
                if qkmer in rkmers:
                    relavent_reads[rid] = set(rkmers)


    def match_lr_kmers(self, relavent_reads):

        for rid, rkmers in self.reads_dict.items():
            for qkmer in self.kmers:
                if qkmer in rkmers:
                    relavent_reads[rid] = set(rkmers)
