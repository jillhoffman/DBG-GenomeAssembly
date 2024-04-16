import unittest
import mapping
import makekmers
import matchkmers
import order
import output
from Bio import SeqIO

rr_overlap = {}
relevant_reads = {}
reads_dict = {'NOMATCH': ['GTGC', 'TGCG', 'GCGA', 'CGAG', 'GAGT', 'AGTC', 'GTCG', 'TCGA', 'CGAG', 'GAGT', 'AGTC', 'GTCG'],
              'FORWARD': ['AGTC', 'GTCG', 'TCGT', 'CGTA', 'GTAT'],
              'BACKWARD': ['CCAG', 'CAGG', 'AGGT', 'GGTG', 'GTGC', 'TGCG']
              }

rr = {
    'FORWARD': ['AGTC', 'GTCG', 'TCGT', 'CGTA', 'GTAT'],
    'BACKWARD': ['CCAG', 'CAGG', 'AGGT', 'GGTG', 'GTGC', 'TGCG']
}
query = [s for s in SeqIO.parse('./../example_data/unittestq.fasta', 'fasta')]
reads = [s for s in SeqIO.parse("./../example_data/unittest.fasta", 'fasta')]


class TestKmers(unittest.TestCase):

    def testMakekmers(self):

        for q in query:
            read_str = str(q.seq)
            test = makekmers.MakeKmers(read_str, 4)
            qkmers = test.make_kmers()
            self.assertEqual(qkmers, ['GTGC', 'TGCG', 'GCGA', 'CGAG', 'GAGT', 'AGTC', 'GTCG'])

    def testMatchkmers(self):

        for q in query:
            read_str = str(q.seq)
            testMake = makekmers.MakeKmers(read_str, 4)
            qkmers = testMake.make_kmers()

            testMatch = matchkmers.MatchKmers(qkmers, reads_dict)
            self.assertEqual(testMatch.match_kmers(relevant_reads, rr_overlap, 'F'), rr)

    def testMakeMaps(self):
        km1FMap = {}
        mF = mapping.map(km1FMap)
        km1All = [('GH', 'HI'), ('CD', 'BC')]

        for km1 in km1All:
            startF = km1[0]
            targetF = km1[1]
            mF.makeMaps(startF, targetF)

        self.assertEqual(km1FMap, {'GH': ['HI'], 'CD':['BC']})

    def testOrderForward(self):
        fsD = {}
        kd = order.Order(2)
        newForQuery = 'CDEFGH'
        mapF = {1: {'GH': ['HI']},
                2: {'HI': ['IJ']}
        }
        jIDs = [1, 2]


        for i in range(len(jIDs)):
            ID = jIDs[i]

            prevForQuery = newForQuery
            newForQuery = kd.forwardSort(prevForQuery, mapF, ID, i, fsD)

        self.assertEqual(newForQuery, 'CDEFGHIJ')

    def testOrderBackward(self):
        bsD = {}
        kd = order.Order(2)
        newForQuery = 'CDEFGH'

        mapB = {1: {'CD': ['BC']},
                2: {'BC': ['AB']}
                }

        jIDs = [1, 2]

        for i in range(len(jIDs)):
            ID = jIDs[i]

            prevForQuery = newForQuery
            newForQuery = kd.backwardSort(prevForQuery, mapB, ID, i, bsD)

        self.assertEqual(newForQuery, 'ABCDEFGH')

    def testOutput(self):

        o = output.outputNeeds(reads)
        sid = 'FORWARD'
        cut = 'TCGT'
        sstart, send = o.ssidInfo(sid, cut)

        self.assertEqual((sstart, send), (2, 6))


