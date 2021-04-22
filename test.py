import unittest
from pileup_reader import pileup_reader, preprocess_bases, get_indel_string
from variant_caller import VariantCaller

class TestPreprocess(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(preprocess_bases(''), '')
        
    def test_normal(self):
        self.assertEqual(preprocess_bases(',.$...,,,.,..^+.'),'.............')
        self.assertEqual(preprocess_bases(',.$$$..A..,'),'....A...')
        self.assertEqual(preprocess_bases('$$$^.'),'')
        
        
class TestIndelString(unittest.TestCase):
    def test_normal(self):
        self.assertEqual(get_indel_string('.+1AAT'),'AAT')
        self.assertEqual(get_indel_string('.-1C'),'C')
        self.assertEqual(get_indel_string('.-2G'),'GG')
        self.assertEqual(get_indel_string('.+3CGC'),'CGCCGCCGC')         
        self.assertEqual(get_indel_string('.-12G'),'GGGGGGGGGGGG')

class TestPileupReader(unittest.TestCase):        
    def test_normal(self):
        pileup_lines = []
        pileup_lines.append({'chromosome': '21', 'position': 9483252,
                        'ref_base': 'T', 'read_count': 1,
                        'read_bases': '.', 'qualities': 'B',
                        'A': 0, 'C': 0, 'G': 0, 'T': 1,
                        'insertions': [], 'deletitions': []})
    
        pileup_lines.append({'chromosome': '21', 'position': 9483266,
                        'ref_base': 'T', 'read_count': 3,
                        'read_bases': '..^<.', 'qualities': '@@>',
                        'A': 0, 'C': 0, 'G': 0, 'T': 3,
                        'insertions': [], 'deletitions': []})
    
        pileup_lines.append({'chromosome': '22', 'position': 41616770,
                        'ref_base': 'G', 'read_count': 90,
                        'read_bases': '...,..,....,,............,.,..........-1A,,,,,,...,.,...,..,.,,..,...........,..,,.,........^],',
                        'qualities': 'fA]@@ZA_HhGGG<rgAJpIkGCoKJIBV=AFIJJJHk@8B@BFIGF>J@eJA9CJ0II>FFBJJ<IBHJJHJJEDH@8AEHJFHCG?AC',
                        'A': 0, 'C': 0, 'G': 89, 'T': 0,
                        'insertions': [], 'deletitions': [['A', 1]]})
    
        i = 0
        for item in pileup_reader('test_data/test.pileup'):
            self.assertEqual(item, pileup_lines[i])
            i += 1
    
        
class TestVariantCaller(unittest.TestCase):
    def test_normal(self):
        variant_caller = VariantCaller()
        mockPositionInfo = { 'A' : 8, 'G' : 1, 'C' : 1, 'T' : 1 , 'ref_base' : 'A'}
        variant_caller.call_variant(mockPositionInfo, 0.8)
        self.assertEqual(mockPositionInfo['alts'],'.')
    
        mockPositionInfo = { 'A' : 1, 'G' : 7, 'C' : 1, 'T' : 8 , 'ref_base' : 'G'}
        variant_caller.call_variant(mockPositionInfo, 0.8)
        self.assertEqual(mockPositionInfo['genotype'],(0, 1))
        self.assertEqual(mockPositionInfo['alts'],['T'])
    
        mockPositionInfo = { 'A' : 1, 'G' : 7, 'C' : 1, 'T' : 8 , 'ref_base' : 'A'}
        variant_caller.call_variant(mockPositionInfo, 0.8)
        self.assertEqual(mockPositionInfo['genotype'],(1, 2))
        self.assertEqual(mockPositionInfo['alts'],['T', 'G'])
    
    def test_indels(self):
        def test_one_insert(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 7, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'G', 'insertions': [('ACAC', 8)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (0, 1))
            self.assertEqual(mockPositionInfo['alts'], ['GACAC'])
            self.assertEqual(mockPositionInfo['ref_base'], 'G')

        def test_one_snv_one_insert(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 7, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'A', 'insertions': [('ACAC', 8)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (1, 2))
            self.assertEqual(mockPositionInfo['alts'], ['AACAC', 'G'])
            self.assertEqual(mockPositionInfo['ref_base'], 'A')
        
        def test_two_inserts(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 1, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'A', 'insertions': [('ACAC', 8), ('GTGT', 7)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (1, 2))
            self.assertEqual(mockPositionInfo['alts'], ['AACAC', 'AGTGT'])
            self.assertEqual(mockPositionInfo['ref_base'], 'A')
        
        def test_one_delete(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 7, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'G', 'deletitions': [('ACAC', 8)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (0, 1))
            self.assertEqual(mockPositionInfo['alts'], ['G'])
            self.assertEqual(mockPositionInfo['ref_base'], 'GACAC')
        
        def test_one_delete_one_snv(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 7, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'A', 'deletitions': [('ACAC', 8)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (1, 2))
            self.assertEqual(mockPositionInfo['alts'], ['A', 'GACAC'])
            self.assertEqual(mockPositionInfo['ref_base'], 'AACAC')
        
        def test_one_delete_one_insert(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 1, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'T', 'deletitions': [('ACAC', 8)], 'insertions' : [('GT', 7)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (1, 2))
            self.assertEqual(mockPositionInfo['alts'], ['T', 'TGTACAC'])
            self.assertEqual(mockPositionInfo['ref_base'], 'TACAC')

        def test_two_deletes(variant_caller):
            mockPositionInfo = { 'A' : 1, 'G': 1, 'C' : 1, 'T' : 1, \
                                'ref_base' : 'T', 'deletitions': [('ACAC', 8), ('AC', 7)]}
            variant_caller.call_variant(mockPositionInfo)
            self.assertEqual(mockPositionInfo['genotype'], (1, 2))
            self.assertEqual(mockPositionInfo['alts'], ['T', 'TAC'])
            self.assertEqual(mockPositionInfo['ref_base'], 'TACAC')
        
        variant_caller = VariantCaller()
        test_one_snv_one_insert(variant_caller)
        test_one_insert(variant_caller)
        test_two_inserts(variant_caller)
        test_one_delete(variant_caller)
        test_one_delete_one_snv(variant_caller)
        test_one_delete_one_insert(variant_caller)
        test_two_deletes(variant_caller)

def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestPreprocess('test_empty'))
    suite.addTest(TestPreprocess('test_normal'))
    suite.addTest(TestIndelString('test_normal'))
    suite.addTest(TestPileupReader('test_normal'))
    suite.addTest(TestVariantCaller('test_normal'))
    suite.addTest(TestVariantCaller('test_indels'))
    return suite

def main():
    runner = unittest.TextTestRunner()
    runner.run(suite())
        
if __name__ == '__main__':
    main()