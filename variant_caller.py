import numpy as np

class VariantCaller(object):
    def __init__(self):
        pass

    def __choose_variant__(self, candidate_base_count, error_probability):
        """
        Given two most common bases b and b', we want to calculate the likelihood of each possible variant
        P(variant | reads) = P(reads | variant) * P(variant) / P(reads)
            => P(variant | reads) is proportional to P(reads | variant) if we ignore P(variant)
        
        If we have k b' bases and n-k b bases, and the error probability for one read base of p:   
            P(reads | variant is b) = C_nk * p^k * (1 - p)^(n-k)
            P(reads | variant is b') = C_nk * p^(n - k) * (1 - p)^(n - k)
            P(reads | variant is bb') = C_nk / 2^n 

        As we have C_nk at all positions, we should just compare the rest. We can do that easier if we take the natural log:
            ln(P(reads | variant is b)) is proportional to k * ln(p) + (n-k) * ln(1-p)
            ln(P(reads | variant is b')) is proportional to (n-k) * ln(p) + k * ln(1-p)
            ln(P(reads | variant is bb')) is proportional to n * ln(1/2)
        """
        first_candidate_base, first_candidate_base_count = candidate_base_count[0]
        second_candidate_base, second_candidate_base_count = candidate_base_count[1]
        n = first_candidate_base_count + second_candidate_base_count
        k = second_candidate_base_count

        first_base_log_likelihood = k * np.log(error_probability) + (n - k) * np.log(1 - error_probability)
        second_base_log_likelihood = (n - k) * np.log(error_probability) + k * np.log(1 - error_probability)
        diploidy_log_likelihood = n * np.log(0.5)
        total_likelihood = np.exp(first_base_log_likelihood) + np.exp(second_base_log_likelihood) + np.exp(diploidy_log_likelihood)

        if first_base_log_likelihood >= second_base_log_likelihood and first_base_log_likelihood >= diploidy_log_likelihood:
            return (first_candidate_base, np.exp(first_base_log_likelihood) / total_likelihood)
        elif second_base_log_likelihood >= first_base_log_likelihood and second_base_log_likelihood >= diploidy_log_likelihood:
            return (second_candidate_base, np.exp(second_base_log_likelihood) / total_likelihood)
        else:
            return (first_candidate_base + second_candidate_base, np.exp(diploidy_log_likelihood) / total_likelihood)
        

    def call_variant(self, genomePositionInfo, error_probability = None):
        if error_probability == None:
            # Deduce error probability from quality
            error_probability = 0.001
        
        base_count = { base: genomePositionInfo[base] for base in {'A', 'G', 'C', 'T'} }
        
        # Only keep two most likely bases
        candidate_bases = sorted(base_count, key=base_count.get, reverse=True)[:2]
        candidate_base_count = [ (base, base_count[base]) for base in candidate_bases ]

        return self.__choose_variant__(candidate_base_count, error_probability)

def main():
    mockPositionInfo = { 'A' : 1, 'G' : 2, 'C' : 1, 'T' : 8 }
    variant_caller = VariantCaller()
    variant, probability = variant_caller.call_variant(mockPositionInfo)
    print(variant, probability)
    assert(variant == 'TG')

    mockPositionInfo = { 'A' : 8, 'G' : 1, 'C' : 1, 'T' : 1 }
    variant, probability = variant_caller.call_variant(mockPositionInfo, 0.2)
    print(variant, probability)
    assert(variant == 'A')

if __name__ == "__main__":
    main()