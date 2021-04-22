import numpy as np
from scipy.stats import binom

class VariantCaller(object):
    def __init__(self):
        pass

    def __calculate_most_probable_variant__(self, candidate_variant_count, correct_probability):
        """
        Parameters
        ----------
        candidate_variant_count: List of (str, int) tuples
            List with two most common variants and their counts
        correct_probability: float
            Probability that one nucleotide in a read is correct
        
        Returns
        -------
        (List, float)
        Chosen variants and estimated correctness probability

        Given two most common variants v and v', we want to calculate the likelihood of each possible variant
        P(variant | reads) = P(reads | variant) * P(variant) / P(reads)
            => P(variant | reads) is proportional to P(reads | variant) if we ignore P(variant)
        
        If we have k v bases and n-k v' bases, and the correct probability for one read position of p:   
            P(reads | variant is v) = C_nk * p^k * (1 - p)^(n-k)
            P(reads | variant is v') = C_nk * p^(n - k) * (1 - p)^k
            P(reads | variant is vv') = C_nk * (1/2)^n

        We can cross out the C_nk as it doesn't affect the comparison. 
        """
        first_candidate_variant, first_candidate_variant_count = candidate_variant_count[0]
        second_candidate_variant, second_candidate_variant_count = candidate_variant_count[1]
        n = first_candidate_variant_count + second_candidate_variant_count
        k = first_candidate_variant_count

        first_variant_likelihood = correct_probability ** k * (1 - correct_probability) ** (n - k)
        second_variant_likelihood = (1 - correct_probability) ** k * correct_probability ** (n - k)
        diploidy_likelihood = (1/2) ** n
        total_likelihood = first_variant_likelihood + second_variant_likelihood + diploidy_likelihood

        if first_variant_likelihood >= second_variant_likelihood and first_variant_likelihood >= diploidy_likelihood:
            return ([first_candidate_variant], first_variant_likelihood / total_likelihood)
        elif second_variant_likelihood >= first_variant_likelihood and second_variant_likelihood >= diploidy_likelihood:
            return ([second_candidate_variant], second_variant_likelihood / total_likelihood)
        else:
            return ([first_candidate_variant, second_candidate_variant], diploidy_likelihood / total_likelihood)
        

    def call_variant(self, genomePositionInfo, correct_probability = 0.8, use_read_quality = False):
        """ Chooses variant and updates genomePositionInfo dictionary with chosen variant genotype, alts field, more

        Parameters
        ----------
        genomePositionInfo: dictionary
            All info about one pileup position
        correct_probability: float
            Probability that one nucleotide in a read is correct
        use_read_quality: bool
            Whether to use correct_probability estimate calculated using read qualities
        """
        if use_read_quality:
            correct_probability = genomePositionInfo['average_quality']
        variant_count = { (base, 'SNV'): genomePositionInfo[base] for base in {'A', 'G', 'C', 'T'} }
        
        # Treat insertions and deletitions the same as SNVs
        if 'insertions' in genomePositionInfo:
            for insertion_string, insertion_count in genomePositionInfo['insertions']:
                variant_count[(insertion_string, 'INS')] = insertion_count
        if 'deletitions' in genomePositionInfo:
            for deletition_string, deletition_count in genomePositionInfo['deletitions']:
                variant_count[(deletition_string, 'DEL')] = deletition_count
        
        # Only keep two most likely bases
        
        candidate_variants = sorted(variant_count, key=variant_count.get, reverse=True)[:2]
        candidate_variant_count = [(variant, variant_count[variant]) for variant in candidate_variants]

        # Check if no candidate variant exists:
        if candidate_variant_count[0][1] == 0:
            genomePositionInfo['vaf'] = 1
            genomePositionInfo['genotype'] = (0, 0)
            genomePositionInfo['alts'] = '.'
            return

        # Check if only one option exists and skip calculations
        if candidate_variant_count[1][1] == 0:
            most_probable_variant = [candidate_variant_count[0][0]]
            confidence = 1
        else:
            most_probable_variant, confidence = self.__calculate_most_probable_variant__(candidate_variant_count, correct_probability)

        ref_variant_present = len([variant[0] for variant in most_probable_variant if (variant[0] == genomePositionInfo['ref_base'] and variant[1] == 'SNV')]) > 0

        alt_variants = [variant[0] for variant in most_probable_variant if variant[0] != genomePositionInfo['ref_base']]

        """
        After calling the variants, we need to store them in the vcf format.
        That means determining the ref_base field, alts field and genotype.

        |------------------------------------------------------------|
        |                    Table for some SNVs                     |
        |------------------------------------------------------------|
        |original ref|called variants|ref in vcf|alts in vcf|genotype|
        |------------------------------------------------------------|
        |     G      |       G       |     G    |     .     |   0/0  |
        |     G      |       T       |     G    |     T     |   1/1  |
        |     G      |      G,T      |     G    |     T     |   0/1  |
        |     G      |      A,T      |     G    |    A,T    |   1/2  |
        --------------------------------------------------------------
        """

        genomePositionInfo['vaf'] = confidence
        if len(alt_variants) == 0:
            genomePositionInfo['genotype'] = (0, 0)
            genomePositionInfo['alts'] = '.'
            return
        if ref_variant_present:
            genomePositionInfo['genotype'] = (0, 1)
        else:
            if len(alt_variants) == 1:
                genomePositionInfo['genotype'] = (1, 1)
            else:
                genomePositionInfo['genotype'] = (1, 2)

        # Extract variant types to diversify INDELs from SNVs
        alt_variant_types = [variant[1] for variant in most_probable_variant if variant[0] != genomePositionInfo['ref_base']]

        """
        Insertions and deletitions require some extra work when writing to alt and ref_base field, 
        given that we might want to change ref_base if we have deletitions, or store insertions as
        ref_base + insertion. Some examples of INDEL storing in vcf:

        |------------------------------------------------------------|
        |                    Table for some INDELS                   |
        |------------------------------------------------------------|
        |original ref|called variants|ref in vcf|alts in vcf|genotype|
        |------------------------------------------------------------|
        |     G      |     CTC IN    |     G    |    GCTC   |   1/1  |
        |     G      |     CTC DEL   |   GCTC   |     G     |   1/1  |
        |     G      |AC DEL, ACA DEL|   GACA   |    GA,G   |   1/2  |
        |     G      | TC DEL, AC IN |   GTC    |  G,GTCAC  |   1/2  |
        --------------------------------------------------------------
        """
        for position in range(len(alt_variants)):
            if alt_variant_types[position] == 'INS':
                alt_variants[position] = genomePositionInfo['ref_base'] + alt_variants[position]

        longest_deletition_string = ''
        for position in range(len(alt_variants)):
            if alt_variant_types[position] == 'DEL':
                if len(alt_variants[position]) > len(longest_deletition_string):
                    longest_deletition_string = alt_variants[position]

        for position in range(len(alt_variants)):
            if alt_variant_types[position] != 'DEL':
                alt_variants[position] += longest_deletition_string
            else:
                shorter_deletition_len = len(alt_variants[position])
                alt_variants[position] = genomePositionInfo['ref_base'] + longest_deletition_string[shorter_deletition_len:]
        genomePositionInfo['ref_base'] = genomePositionInfo['ref_base'] + longest_deletition_string

        genomePositionInfo['alts'] = alt_variants
        

def main():
    variant_caller = VariantCaller()

    mockPositionInfo = { 'A' : 8, 'G' : 1, 'C' : 1, 'T' : 1 , 'ref_base' : 'A'}
    variant_caller.call_variant(mockPositionInfo, 0.2)
    assert(mockPositionInfo['alt'] == '.')

    mockPositionInfo = { 'A' : 1, 'G' : 2, 'C' : 1, 'T' : 8 , 'ref_base' : 'G'}
    variant_caller.call_variant(mockPositionInfo)
    assert(mockPositionInfo['genotype'] == (0, 1))
    assert(mockPositionInfo['alt'] == ['T'])

    mockPositionInfo = { 'A' : 1, 'G' : 2, 'C' : 1, 'T' : 8 , 'ref_base' : 'A'}
    variant_caller.call_variant(mockPositionInfo)
    assert(mockPositionInfo['genotype'] == (1, 2))
    assert(mockPositionInfo['alt'] == ['T', 'G'])


if __name__ == "__main__":
    main()