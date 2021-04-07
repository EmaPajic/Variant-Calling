from pileup_reader import pileup_reader
from variant_caller import VariantCaller
from vcf_writer import create_vcf_file, write_vcf_line
import time
import argparse

def main():
    parser = argparse.ArgumentParser(description='Runs variant calling on pileup file and stores in vfc file')
    parser.add_argument('--use-read-quality', default='False', action='store_true',
                        help='tells the algorithm to estimate p from read qualities')
    parser.add_argument('--call-less-positions', default=False, action='store_true',
                        help='tells the program to call less positions (not whole pileup file)')
    parser.add_argument('--input-file', default='merged-normal.pileup', type=str,
                        help='path to input file in pileup format')
    parser.add_argument('--output-file', default='Make name from input name', type=str,
                        help='name for the output vcf file. If not given, will be created from input file name')
    parser.add_argument('--p', default='0.99', type=float,
                        help='probability estimate of one nucleotide read being correct, used by vc algorithm')
    parser.add_argument('--positions-to-call', default='10000', type=int,
                        help='how many positions to call if call-less-positions set to true')
    args = parser.parse_args()
    if args.output_file == 'Make name from input name':
        args.output_file = args.input_file + '.vcf'
    
    variant_caller = VariantCaller()
    sample = 'SAMPLE1'

    create_vcf_start = time.time()
    vcf = create_vcf_file(args.output_file, sample)
    create_vcf_end = time.time()
    print('Vcf header created. Elapsed time: {}'.format(create_vcf_end - create_vcf_start))

    main_loop_start = time.time()
    position_count = 0
    variant_caller_time = 0
    positions_with_variants = 0
    write_vcf_time = 0

    for pileup_line in pileup_reader(args.input_file):
        variant_caller_start = time.time()
        variant_caller.call_variant(pileup_line, args.p, args.use_read_quality)
        if pileup_line['alts'] != '.':
            positions_with_variants += 1
        variant_caller_time += time.time() - variant_caller_start

        write_vcf_start = time.time()
        write_vcf_line(pileup_line, vcf, sample)
        write_vcf_time = time.time() - write_vcf_start

        position_count += 1
        if args.call_less_positions and (position_count >= args.positions_to_call):
            break
    
    main_loop_end = time.time()
    total_running_time = main_loop_end - main_loop_start

    print('Processed {} positions. Found variants at {} positions.'.format(position_count, positions_with_variants))

    print('Total running time is {}'.format(total_running_time))
    print('Pileup reader: {}'.format(total_running_time - variant_caller_time - write_vcf_time))
    print('Variant calling: {}'.format(variant_caller_time))
    print('Vcf writing: {}'.format(write_vcf_time))

if __name__ == '__main__':
    main()