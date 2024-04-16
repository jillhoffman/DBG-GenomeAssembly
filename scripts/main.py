from Bio import SeqIO
import makekmers
import matchkmers
import order
import output
import time
import mapping

def main():

    # make_kmer_st = time.time()
    # reads = [s for s in SeqIO.parse("./../example_data/READS.fasta", 'fasta')]
    # initial_queries = [r for r in SeqIO.parse("./../example_data/QUERY.fasta", 'fasta')]
    # kmer_size = 4

    fasta_queries = sys.argv[1]
    initial_queries = [r for r in SeqIO.parse(fasta_queries,'fasta')]
    fasta_reads = sys.argv[2]
    reads = [s for s in SeqIO.parse(fasta_reads,'fasta')]
    kmer_size = int(sys.argv[3])

    reads_dict = {}


    print(f"creating {kmer_size}-mers...")

    for read in reads:
        read_str = str(read.seq)
        rkmercuts = makekmers.MakeKmers(read_str, kmer_size)
        initial_r_kmers = rkmercuts.make_kmers()
        reads_dict[read.id] = initial_r_kmers

    num_reads = len(reads_dict)

    for query in initial_queries:
        qid = query.id
        createdContigs = {}
        query_str = str(query.seq)
        qkmercuts = makekmers.MakeKmers(query_str, kmer_size)
        initial_q_kmers = qkmercuts.make_kmers()

        print(f"matching {num_reads} reads to initial query...")

        relavent_reads = {}
        rr_overlap = {}
        km1_size = kmer_size - 1
        dbgForwardMap = {}
        dbgBackwardMap = {}

        matching = matchkmers.MatchKmers(initial_q_kmers, reads_dict)
        matching.match_kmers(relavent_reads, rr_overlap, "T")
        query_rr = len(relavent_reads)
        jIDs = []

        max_overlap = max(rr_overlap.values())
        id_max_overlap = list(filter(lambda x: rr_overlap[x] == max_overlap, rr_overlap))
        submax = id_max_overlap[:5]

        for id in submax:
            rrkmers = relavent_reads[id]
            rrmatching = matchkmers.MatchKmers(rrkmers, reads_dict)
            rrmatching.match_kmers(relavent_reads, rr_overlap, "F")

        if len(rr_overlap.values()) != 0:

            for rrID in dict(sorted(rr_overlap.items(), key=lambda x:x[1], reverse=True)):

                km1All = []
                jIDs.append(rrID)
                rrKmers = relavent_reads[rrID]

                for rrKmer in rrKmers:
                    km1cuts = makekmers.MakeKmers(rrKmer, km1_size)
                    km1s = km1cuts.make_km1ers(initial_q_kmers)

                    if len(km1s) > 0:
                        km1All.append(km1s)

                km1FMap = {}
                km1BMap = {}
                mF = mapping.map(km1FMap)
                mB = mapping.map(km1BMap)

                for km1 in km1All:

                    startF = km1[0]
                    targetF = km1[1]

                    mF.makeMaps(startF, targetF)
                    mB.makeMaps(targetF, startF)

                dbgForwardMap[rrID] = km1FMap
                dbgBackwardMap[rrID] = km1BMap

            print(f"matched {query_rr} reads")

        else:
            print(f"no reads were matched, try a smaller k-mer size")

        print("time elapsed: %s seconds" % (round(time.time() - make_kmer_st)))
        print('extending original query...')
        extSt = time.time()
        kd = order.Order(km1_size)

        for i in range(100):
            md = {}
            c, fsD, bsD = DBG(query_str, kd, dbgForwardMap, dbgBackwardMap, jIDs)
            md['length'] = len(c)
            md['contig']=c
            md['forward']=fsD
            md['backward']=bsD

            createdContigs[i] = md

        maxL = max(val['length'] for val in createdContigs.values())

        longestMatches = [c for id, c in createdContigs.items() if c['length'] == maxL]
        numLongMatches = len(longestMatches)

        print(f"assembled {numLongMatches} alleles with {maxL} bps")

        alleles = open("../output/ALLELES.fasta", "w")

        alnData = []
        for i,c in enumerate(longestMatches):

            o = output.outputNeeds(reads)

            cid = f'CONTIG{i}_L={c["length"]}'
            alleles.write(">" + cid + "\n" + c["contig"] + "\n")

            backSort = c['backward']
            forSort = c['forward']
            contig = c['contig']
            lenB = len(backSort)
            totLen = 0

            if lenB == 1:
                sid = list(backSort.keys())[0]
                qstart = 0
                qend = len(backSort[sid])
                cut = contig[qstart:qend]

                sstart, send = o.ssidInfo(sid, cut)

                alnLine = f'{sid}\t{cid}\t{sstart+1}\t{send}\t{qstart+1}\t{qend}\n'
                alnData.append(alnLine)



                totLen += qend

            elif lenB > 1:

                for i in range(lenB + 1):

                    next = lenB - i
                    nextID = list(backSort.keys())[next]

                    lenNext = len(backSort[nextID])

                    qstart = totLen
                    qend = totLen + lenNext

                    cut = contig[qstart:qend]
                    sstart, send = o.ssidInfo(nextID, cut)

                    alnLine = f'{nextID}\t{cid}\t{sstart + 1}\t{send}\t{qstart + 1}\t{qend}\n'
                    alnData.append(alnLine)

                    totLen += lenNext

            sstart = 1
            send = len(query_str)
            qstart = totLen
            qend = totLen + len(query_str)
            totLen += len(query_str)

            alnLine = f'{qid}\t{cid}\t{sstart}\t{send}\t{qstart + 1}\t{qend}\n'
            alnData.append(alnLine)

            if len(forSort) == 1:
                sid = list(forSort.keys())[0]
                qstart = totLen
                qend = totLen + len(forSort[sid])
                cut = contig[qstart:qend]

                sstart, send = o.ssidInfo(sid, cut)

                alnLine = f'{sid}\t{cid}\t{sstart + 1}\t{send}\t{qstart + 1}\t{qend}\n'
                alnData.append(alnLine)


            elif len(forSort) > 1:
                for ID, pos in forSort.items():
                    length = len(pos)
                    qstart = totLen
                    qend = totLen + length

                    alnLine = f'{sid}\t{cid}\t{sstart + 1}\t{send}\t{qstart + 1}\t{qend}\n'
                    alnData.append(alnLine)

                    totLen += length + 1

        alleles.close()
        alnAlleles = open("../output/ALLELES.aln", "w")
        alnAlleles.write("sseqid\tqsid\tsstart\tsend\tqstart\tqend\n")

        for m in alnData:

            alnAlleles.write(m)

        alnAlleles.close()
        print("time elapsed: %s seconds" % (round(time.time() - extSt)))



def DBG(query_str, kd, dbgForwardMap, dbgBackwardMap, jIDs):

    fsD = {}
    bsD = {}

    i=0
    x=0
    newForQuery = query_str
    while x < len(jIDs):
        ID = jIDs[x]

        if i==0:
            prevForQuery = newForQuery
            newForQuery = kd.forwardSort(prevForQuery, dbgForwardMap, ID, i, fsD)

            if len(newForQuery) - len(prevForQuery) == 1:
                i+=1

            else:
                x+=1

        elif len(newForQuery) - len(prevForQuery) == 1:

            prevForQuery = newForQuery
            newForQuery = kd.forwardSort(prevForQuery, dbgForwardMap, ID, i, fsD)
            i+=1

        elif len(newForQuery) - len(prevForQuery) == 0:
            x+=1

    i = 1
    x = 0

    while x < len(jIDs):
        ID = jIDs[x]

        if i==1:
            prevBackQuery = newForQuery
            newBackQuery = kd.backwardSort(prevBackQuery, dbgBackwardMap, ID, i, bsD)
            i += 1

        elif len(newBackQuery) - len(prevBackQuery) == 1:
            prevBackQuery = newBackQuery
            newBackQuery = kd.forwardSort(newBackQuery, dbgBackwardMap, ID, i, fsD)
            i += 1

        elif len(newBackQuery) - len(prevBackQuery) == 0:
            x += 1

    return newBackQuery, fsD, bsD

if __name__ == '__main__':
    main()






