import numpy as np

class VariantCaller(object):
    def __init__(self):
        pass

    def __calculate_most_probable_variant__(self, candidate_variant_count, error_probability):
        """
        Given two most common variants v and v', we want to calculate the likelihood of each possible variant
        P(variant | reads) = P(reads | variant) * P(variant) / P(reads)
            => P(variant | reads) is proportional to P(reads | variant) if we ignore P(variant)
        
        If we have k v' bases and n-k v bases, and the error probability for one read position of p:   
            P(reads | variant is v) = C_nk * p^k * (1 - p)^(n-k)
            P(reads | variant is v') = C_nk * p^(n - k) * (1 - p)^(n - k)
            P(reads | variant is vv') = C_nk / 2^n 

        As we have C_nk at all positions, we should just compare the rest. We can do that easier if we take the natural log:
            ln(P(reads | variant is v)) is proportional to k * ln(p) + (n-k) * ln(1-p)
            ln(P(reads | variant is v')) is proportional to (n-k) * ln(p) + k * ln(1-p)
            ln(P(reads | variant is vv')) is proportional to n * ln(1/2)
        """
        first_candidate_variant, first_candidate_variant_count = candidate_variant_count[0]
        second_candidate_variant, second_candidate_variant_count = candidate_variant_count[1]
        n = first_candidate_variant_count + second_candidate_variant_count
        k = second_candidate_variant_count

        first_variant_log_likelihood = k * np.log(error_probability) + (n - k) * np.log(1 - error_probability)
        second_variant_log_likelihood = (n - k) * np.log(error_probability) + k * np.log(1 - error_probability)
        diploidy_log_likelihood = n * np.log(0.5)
        total_likelihood = np.exp(first_variant_log_likelihood) + np.exp(second_variant_log_likelihood) + np.exp(diploidy_log_likelihood)

        if first_variant_log_likelihood >= second_variant_log_likelihood and first_variant_log_likelihood >= diploidy_log_likelihood:
            return ([first_candidate_variant], np.exp(first_variant_log_likelihood) / total_likelihood)
        elif second_variant_log_likelihood >= first_variant_log_likelihood and second_variant_log_likelihood >= diploidy_log_likelihood:
            return ([second_candidate_variant], np.exp(second_variant_log_likelihood) / total_likelihood)
        else:
            return ([first_candidate_variant, second_candidate_variant], np.exp(diploidy_log_likelihood) / total_likelihood)
        

    def call_variant(self, genomePositionInfo, error_probability = None):
        if error_probability == None:
            # Deduce error probability from quality
            error_probability = 0.001
        
        variant_count = { base: genomePositionInfo[base] for base in {'A', 'G', 'C', 'T'} }
        if 'insertions' in genomePositionInfo:
            for insertion_string, insertion_count in genomePositionInfo['insertions']:
                variant_count[insertion_string] = insertion_count
        if 'deletitions' in genomePositionInfo:
            for deletition_string, deletition_count in genomePositionInfo['deletitions']:
                variant_count[deletition_string] = deletition_count
        
        # Only keep two most likely bases
        candidate_variants = sorted(variant_count, key=variant_count.get, reverse=True)[:2]
        candidate_variant_count = [(variant, variant_count[variant]) for variant in candidate_variants]

        most_probable_variant, confidence = self.__calculate_most_probable_variant__(candidate_variant_count, error_probability)
        
        ref_base_present = len([base for base in most_probable_variant if base == genomePositionInfo['ref_base']]) == 1
        alt_bases = [base for base in most_probable_variant if base != genomePositionInfo['ref_base']]
        
        genomePositionInfo['vaf'] = confidence
        if len(alt_bases) == 0:
            genomePositionInfo['genotype'] = (0, 0)
            genomePositionInfo['alts'] = '.'
            return
        if ref_base_present:
            genomePositionInfo['genotype'] = (0, 1)
        else:
            if len(alt_bases) == 1:
                genomePositionInfo['genotype'] = (1, 1)
            else:
                genomePositionInfo['genotype'] = (1, 2)

        genomePositionInfo['alts'] = alt_bases
        

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