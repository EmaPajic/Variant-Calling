from pileup_reader import pileup_reader
from variant_caller import VariantCaller
from vcf_writer import create_vcf_file, write_vcf_line
import time

def main():
    variant_caller = VariantCaller()
    sample = 'SAMPLE1'
    vcf = create_vcf_file('merged-normal.vcf', sample)
    
    start = time.time()

    for pileup_line in pileup_reader('merged-normal.pileup'):
        variant_caller.call_variant(pileup_line)
        write_vcf_line(pileup_line, vcf, sample)
    
    print(time.time() - start)
if __name__ == '__main__':
    main()