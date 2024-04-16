class outputNeeds:

    def __init__(self, reads):
        self.reads = reads

    def ssidInfo(self, sid, cut):

        for read in self.reads:
            if read.id == sid:
                readSeq = read.seq

                sstart = readSeq.find(cut)
                send = sstart + len(cut)

                return sstart, send
