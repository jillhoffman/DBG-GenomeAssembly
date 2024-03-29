import unittest
import makekmers
import matchkmers
from Bio import SeqIO

relevant_reads = {}
class TestKmers(unittest.TestCase):

    def testMakekmers(self):

        reads = [s for s in SeqIO.parse('./../example_data/unittest.fasta', 'fasta')]
        for read in reads:
            read_str = str(read.seq)
            test = makekmers.MakeKmers(read_str, 3)
            self.assertEqual(test.make_kmers().sort(), ['GGG', 'GGC', 'GCC', 'CCC', 'CCC', 'CCA', 'CAG', 'AGC', 'GCA', 'CAC', \
                                                 'ACC', 'CCT', 'CTG', 'TGA', 'GAC', 'ACT', 'CTT', 'TTG', 'TGG', 'GGC', 'GCC', 'CCT', 'CTC', \
                                                 'TCC', 'CCT', 'CTC', 'TCC', 'CCC'].sort())

    def testMatchkmers(self):
        query = [s for s in SeqIO.parse('./../example_data/QUERY.fasta', 'fasta')]
        reads_dict = {'2S43D:05027:10155': ['GGG', 'GGC', 'GCC', 'CCC', 'CCC', 'CCA', 'CAG', 'AGC', 'GCA', 'CAC', \
                                            'ACC', 'CCT', 'CTG', 'TGA', 'GAC', 'ACT', 'CTT', 'TTG', 'TGG', 'GGC', 'GCC', 'CCT', 'CTC', \
                                            'TCC', 'CCT', 'CTC', 'TCC', 'CCC'],
                      '2S43D:03036:10666': ['GAA', 'AAC', 'ACG', 'CGA', 'GAA', 'AAC', 'ACA', 'CAG', 'AGT', 'GTC', 'TCA', 'CAC', 'ACA', 'CAG', 'AGG', 'GGC', \
                                            'GCA', 'CAA', 'AAA', 'AAT', 'ATC', 'TCC', 'CCA', 'CAG', 'AGG', 'GGC', 'GCG', 'CGG', 'GGT', 'GTC', 'TCA', 'CAT', 'ATG', \
                                            'TGG', 'GGC', 'GCC', 'CCT', 'CTC', 'TCG', 'CGA', 'GAC', 'ACT', 'CTG', 'TGG', 'GGG', 'GGC', 'GCT', 'CTG', 'TGG', 'GGG', 'GGG', \
                                            'GGG', 'GGG', 'GGG', 'GGG', 'GGA', 'GAT', 'ATT', 'TTG', 'TGA', 'GAA', 'AAG', 'AGC', 'GCT', 'CTC', 'TCA', 'CAA', 'AAA', 'AAG', \
                                            'AGT', 'GTA', 'TAT', 'ATA', 'TAA', 'AAA', 'AAA', 'AAG', 'AGC', 'GCG', 'CGT', 'GTC', 'TCG', 'CGC', 'GCC', 'CCC', 'CCG', \
                                            'CGC', 'GCT', 'CTG', 'TGG', 'GGA', 'GAA', 'AAG', 'AGG', 'GGC', 'GCT', 'CTG', 'TGC', 'GCG', 'CGG', 'GGC', 'GCG', 'CGC', \
                                            'GCC', 'CCG', 'CGT', 'GTG', 'TGG', 'GGG', 'GGT', 'GTC', 'TCC', 'CCA', 'CAG', 'AGA', 'GAG', 'AGG', 'GGC', 'GCC', 'CCT', 'CTG', 'TGC', 'GCG', 'CGA', 'GAG', 'AGC', 'GCC', 'CCC', 'CCA', 'CAG', 'AGG', 'GGG', 'GGA'],
                      }

        for q in query:
            qStr = str(q.seq)
            qkmersob = makekmers.MakeKmers(qStr, 3)
            qkmers = qkmersob.make_kmers()

            test = matchkmers.MatchKmers(qkmers, reads_dict)

            self.assertEqual(list(test.match_lr_kmers(relevant_reads).keys())[0], '2S43D:03036:10666')
