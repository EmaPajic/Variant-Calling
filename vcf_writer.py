import pysam
import datetime

def create_vcf_file(path, sample):
    """ Creates VCF header and Variant File. 
    Writes VCF header in Variant File and returns it. 
        
    Parameters
    ----------
    path: str
        Name and path of an output vcf file, for example output/out.vcf
    sample: str
        Name of a sample to add to the VCF file
    
    Returns
    -------
    pysam.VariantFile
        Created VCF file with header written in it
    """
    
    vcf_header = pysam.VariantHeader()
    vcf_header.add_sample(sample)
    
    current_time = datetime.datetime.now()
    date = current_time.strftime('%Y%m%d')
    vcf_header.add_line('##fileDate=' + date)
    vcf_header.add_line('##source=Ema&Nikola')                    
    
    faifile = open("test_data/human_g1k_v37_decoy.fasta.fai")
    for line in faifile:
            split_line = line.split("\t")
            contig = '##contig=<ID=' + str(split_line[0]) + ', length=' + str(split_line[1]) + '>'
            vcf_header.add_line(contig)
    faifile.close()
    vcf_header.add_line("##ALT=<ID=*,Description=Different allele than referent.>")
    vcf_header.add_line("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>")
    vcf_header.add_line("##FORMAT=<ID=VAF,Number=1,Type=String,Description=Variant allele frequency>")
 
    vcf = pysam.VariantFile(path, 'w', header = vcf_header)
    return vcf

def write_vcf_line(pileup_record, vcf, sample): 
    """ Writes a line from a pileup record for a given sample to given vcf file.
    
    Parameters
    ----------
    pileup_record: str
        Line of a pileup file to be written to VCF
    vcf: pysam.VariantFile
        VCF file where to write
    sample: str
        Name of a sample in the VCF file
    
    """
    
    record = vcf.header.new_record()
    record.contig = pileup_record['chromosome']
    record.pos = pileup_record['position']
    record.ref = pileup_record['ref_base']
    record.alts = pileup_record['alts']
    record.samples[sample]['GT'] = pileup_record['genotype']
    record.samples[sample]['VAF'] = str(pileup_record['vaf'])
    
    vcf.write(record)

if __name__ == '__main__':
    sample = 'SAMPLE1'
    vcf = create_vcf_file('vcffile.txt', sample)
    vcf.close()