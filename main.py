from pileup_reader import pileup_reader
from variant_caller import VariantCaller

def main():
    variant_caller = VariantCaller()
    cnt = 0
    for pileup_line in pileup_reader('merged-normal.pileup'):
        variant_caller.call_variant(pileup_line)
        
        if pileup_line['alts'] != '.':
            print(pileup_line)
            cnt += 1

        if cnt == 100:
            break

if __name__ == '__main__':
    main()