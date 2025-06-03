import RNA
import unittest

from drtransformer.rnafolding import call_findpath

class PlatformFindpathTests(unittest.TestCase):
    # NOTE: The path returned by findpath is sometimes different on ubuntu and macos!
    def test_findpath(self):
        seq = 'GCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUU'
        ss1 = '(......((((((((......)))).......))))..)'
        ss2 = '((((.......((((......)))).......)))...)'
        w = 4 * RNA.bp_distance(ss1, ss2)
        path, barrier = call_findpath(seq, ss1, ss2, RNA.md(), w)
        assert barrier == 810
        assert path == [('(......((((((((......)))).......))))..)', 400),
                        ('(.......(((((((......)))).......)))...)', 600),
                        ('(.......((.((((......))))........))...)', 860),
                        ('(........(.((((......))))........)....)', 1210),
                        ('(..........((((......)))).............)', 680),
                        ('((.........((((......)))).........)...)', 770),
                        ('(((........((((......))))........))...)', 430),
                        ('((((.......((((......)))).......)))...)', 280)]
        assert barrier == max((en for ss, en in path)) - path[0][1]
        print(f'\n{seq}\n{ss1}\n{ss2}')
        for (ss, en) in path:
            print(f'{ss} {en:>5d}')

    def test_findpath_split_01(self):
        seq = 'AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUACAAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC'
        ss1 = '...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((..((.(((((......)))))..))..)))...))))............'
        ss2 = '...(((((((.......((((......)))).......))))...)))...............(((((((.......((((......)))).......))))...)))............'
        w = 4 * RNA.bp_distance(ss1, ss2)
        path, barrier = call_findpath(seq, ss1, ss2, RNA.md(), w)
        assert barrier == 410
        assert barrier == max((en for ss, en in path)) - path[0][1]
        print(f'\n{seq}\n{ss1}\n{ss2}')
        for (ss, en) in path:
            print(f'{ss} {en:>5d}')

    def test_findpath_split_02(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.(.((...)).)..)))......(((.......)))...."
        w = 4 * RNA.bp_distance(ss1, ss2)
        path, barrier = call_findpath(seq, ss1, ss2, RNA.md(), w)
        assert barrier == 110
        print(f'\n{seq}\n{ss1}\n{ss2}')
        for (ss, en) in path:
            print(f'{ss} {en:>5d}')


