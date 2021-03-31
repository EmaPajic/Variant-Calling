import numpy as np
import re
from collections import Counter

def preprocess_bases(read_bases):
    """ Returns read without irrelevant characters
    """
    read_bases = read_bases.upper().replace(',', '.')
    read_bases = re.sub('\^.', '', read_bases)
    read_bases = re.sub('\$|\*','',read_bases)
    return read_bases

def count_bases(read_bases):
    """ Returns count of every base
    """
    base_counter = Counter(read_bases).most_common() 
    return base_counter

def get_average_quality(qualities):
    """ Calculates average quality as Phred quality score
    """
    sum_quality = 0
    
    for q in qualities:
        sum_quality += 1 - 10 **-((ord(q) - 33) / 10.0)
        
    return sum_quality / len(qualities)

def get_indel_string(read_bases):
    """ Returns actual indel string
    """
    indel_string = read_bases[2:]
    
    ind_num = [ind.end() for ind in (re.finditer(r'[A-Za-z\.,][0-9]*[^A-Za-z]', read_bases))]

    if len(ind_num) > 0:
        indel_string = read_bases[2:ind_num[0]]
        
        if len(ind_num) != 1:
            for i in range(0,len(ind_num)-1):
                indel_string += int(float(read_bases[ind_num[i]])) * \
                                read_bases[ind_num[i]+1:ind_num[i+1]]
        indel_string += int(float(read_bases[ind_num[-1]])) * \
                        read_bases[ind_num[-1] + 1:]

    return indel_string


def pileup_reader(path):
    """ Reads pileup file, removes irrelevant characters from read, 
        counts bases, detects insertions and deletions and returnes 
        a dictionary with all relevant information.
    """
    with open(path, 'r') as pileup_file:
        i = 0
        for line in pileup_file.readlines():
            split_line = line.split()
            
            pileup_line = {}
            pileup_line['chromosome'] = split_line[0]
            pileup_line['position'] = split_line[1]
            pileup_line['ref_base'] = split_line[2]
            pileup_line['read_count'] = int(split_line[3])
            pileup_line['read_bases'] = split_line[4]
            pileup_line['qualities'] = split_line[5]
            
            pileup_line['average_quality'] = get_average_quality(split_line[5])
            
            pileup_line['A'] = 0
            pileup_line['C'] = 0
            pileup_line['G'] = 0
            pileup_line['T'] = 0
            
            read_bases = preprocess_bases(pileup_line['read_bases'])
            
            ins = re.findall(r'[\.][+][ACGT]*[0-9]*[ACGT]*[0-9]*[ACGT]*', read_bases)
            dels = re.findall(r'[\.][-][ACGT]*[0-9]*[ACGT]*[0-9]*[ACGT]*', read_bases)
            var_insertion = []
            var_deletition = []
            insertion_variants = list(set(ins))
            deletition_variants = list(set(dels))
            
            insertion_variants1 = [get_indel_string(var) for var in insertion_variants]
            deletition_variants1 = [get_indel_string(var) for var in deletition_variants]
            
            var_counts_insertion = [ins.count(indel) for indel in insertion_variants]
            for i in range(0, len(insertion_variants1)):
                var_insertion.append([insertion_variants1[i], var_counts_insertion[i]])
                
            var_counts_deletition = [dels.count(indel) for indel in deletition_variants]
            for i in range(0, len(deletition_variants1)):
                var_deletition.append([deletition_variants1[i], var_counts_deletition[i]])
            
            insertion_variants.sort(key = len, reverse = True)
            for s in insertion_variants:
                read_bases = read_bases.replace(s,'')
                
            deletition_variants.sort(key = len, reverse = True)
            for s in deletition_variants:
                read_bases = read_bases.replace(s,'')
                
            pileup_line['insertions'] = var_insertion
            pileup_line['deletitions'] = var_deletition
            
            read_bases = read_bases.replace('.', pileup_line['ref_base'])
            base_counter = count_bases(read_bases)
            
            for base in base_counter:
                pileup_line[base[0]] = base[1]
            
            yield pileup_line
            
                
if __name__ == '__main__':
    for item in pileup_reader('/home/user/Desktop/EmaPajic/Variant-Calling/merged-normal.pileup'):
        print(item)