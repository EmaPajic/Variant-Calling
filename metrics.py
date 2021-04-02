import numpy as np
import pysam

def precision(tp, fp, fn):
    return tp / (tp + fp)

def recall(tp, fp, fn):
    return tp / (tp + fn)

def f1_score(tp, fp, fn):
    prec = precision(tp, fp, fn)
    rec = recall(tp, fp, fn)
    return 2 * prec * rec / (prec + rec)

def accuracy(tp, fp, fn, tn):
    return (tp + tn) / (tp + fp + tn + fn)

def mcc(tp, fp, fn, tn):
    return 1.0
    return (tp * tn - fp * fn) /\
            np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
          
            
def get_statistics(bfctools_vcf_file, vcf_file):        
    bcftools_vcf = pysam.VariantFile(bcftools_vcf_file, "r")
    vcf = pysam.VariantFile(vcf_file, "r")
        
    tp = 0
    fp = 0
    fn = 0
    tn = 0
    
    data_bcftools = {}
    for record in bcftools_vcf.fetch():
        #print(record.samples['HCC1143BL']['GT'])
        if record.samples['HCC1143BL']['GT'] != (0, 0):
            data_bcftools[(record.chrom, record.pos)] = \
            [record.ref, record.alts, record.samples['HCC1143BL']['GT']]
        
    for record in vcf.fetch():
        #print(record.samples['SAMPLE1']['GT'])
        genotype = record.samples['SAMPLE1']['GT']
        
        if genotype == (0, 0) and (record.chrom, record.pos) not in data_bcftools:
            tn += 1
        elif genotype == (0, 0) and (record.chrom, record.pos) not in data_bcftools:
            fn += 1
        elif genotype != (0, 0) and (record.chrom, record.pos) in data_bcftools:
            fp += 1
        else:
            tp += 1
    
    bcftools_vcf.close()
    vcf.close()
    return tp, fp, fn, tn
    
if __name__ == '__main__':
    bcftools_vcf_file = "merged-normal.bam.mpileup.vcf.called.vcf"
    vcf_files = ["merged-normal.vcf"]
    
    TP = []
    FP = []
    FN = []
    TN = []
    precision_list = []
    recall_list = []
    f1_score_list = []
    accuracy_list = []
    mcc_list = []
    
    for vcf_file in vcf_files:
        tp, fp, fn, tn = get_statistics(bcftools_vcf_file, vcf_file)
        TP.append(tp)
        FP.append(fp)
        FN.append(fn)
        TN.append(tn)
        precision_list.append(precision(tp, fp, fn))
        recall_list.append(recall(tp, fp, fn))
        f1_score_list.append(f1_score(tp, fp, fn))
        accuracy_list.append(accuracy(tp, fp, fn, tn))
        mcc_list.append(mcc(tp, fp, fn, tn))
        
    print(precision_list)
    print(recall_list)
    print(f1_score_list)
    print(accuracy_list)
    print(mcc_list)
          
    
    
        
        
        