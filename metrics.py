import numpy as np
import pysam
import matplotlib.pyplot as plt

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
    return (tp * tn - fp * fn) /\
            np.sqrt(1.0 * (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
          
            
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
    vcf_files = ["merged-normal.80.vcf", "merged-normal.99.vcf"]
    
    p = []
    
    for vcf_file in vcf_files:
        info = vcf_file.split('.')
        p.append(float(info[1]) / 100)
    
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
       
    for i in range(len(vcf_files)):
        print('Probability: {}'.format(p[i]))
        print('Precision: {}'.format(precision_list[i]))
        print('Recall: {}'.format(recall_list[i]))
        print('F1 score: {}'.format(f1_score_list[i]))
        print('Accuracy: {}'.format(accuracy_list[i]))
        print('MCC score: {}'.format(mcc_list[i]))
        print('')
    
    plt.figure('Precision')
    plt.title('Precision')
    plt.xlabel('Probability')
    plt.ylabel('Precision')
    plt.plot(p, precision_list)
    
    plt.figure('Recall')
    plt.title('Recall')
    plt.xlabel('Probability')
    plt.ylabel('Recall')
    plt.plot(p, recall_list)
    
    plt.figure('F1 score')
    plt.title('F1 score')
    plt.xlabel('Probability')
    plt.ylabel('F1 score')
    plt.plot(p, f1_score_list)
    
    plt.figure('Accuracy')
    plt.title('Accuracy')
    plt.xlabel('Probability')
    plt.ylabel('Accuracy')
    plt.plot(p, accuracy_list)
    
    plt.figure('MCC score')
    plt.title('MCC score')
    plt.xlabel('Probability')
    plt.ylabel('MCC score')
    plt.plot(p, mcc_list)
    
          
    
    
        
        
        