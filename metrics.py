import numpy as np
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn


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
    """ Calculates the number of true positives, false positives, 
    false negatives and true negatives. True positives represent variants 
    appearing both in bfctools VCF file and ours, false positives are variants
    appearing in our VCF file but not in bfctools, false negatives are 
    variants appearing in bfctools VCF file but not in ours, and finally, 
    true negatives are not appearing in both.
    
    Parameters
    ----------
    bfctools_vcf_file: pysam.VariantFile
        VCF file created by bfctools call tool
    vcf_file: pysam.VariantFile
        VCF file created by our algorithm
        
    Returns
    -------
    (int, int, int, int)
        Number of true positives, false positives, false negatives and
        true negatives
    """    
    
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
        elif genotype == (0, 0) and (record.chrom, record.pos) in data_bcftools:
            fn += 1
        elif genotype != (0, 0) and (record.chrom, record.pos) not in data_bcftools:
            fp += 1
        else:
            tp += 1
    
    bcftools_vcf.close()
    vcf.close()
    return tp, fp, fn, tn
    
def metrics(bcftools_vcf_file, vcf_files):
    """ Prints precision, recall, F1 score, accuracy, MCC score and 
    confusion matrix for each VCF file. Plots precision, recall, F1 score,
    accuracy and MCC score against probabilities. 
    
    Parameters
    ----------
    bfctools_vcf_file: str
        path to VCF file created by bfctools call tool
    vcf_files: list of str
        paths to VCF files created by our algorithm with different 
        probabilities
        
    """  
    
    p = []
    
    for vcf_file in vcf_files:
        info = vcf_file.split('.')
        p.append(float(info[2]) / 100)
    
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
        
        confusion_matrix = np.array([[TN[i], FP[i]],[FN[i], TP[i]]])
        df_cm = pd.DataFrame(confusion_matrix, range(2), range(2))
        plt.figure(i)
        sn.set(font_scale=1.4)
        ax = sn.heatmap(df_cm, annot=True, annot_kws={"size": 16}, fmt="d", cmap="YlGnBu")
        
        ax.set(xlabel='Predicted', ylabel='True')

        plt.show() 
    
    plt.figure('Metrics')
    plt.title('Metrics')
    plt.xlabel('Probability')
    plt.ylabel('Metrics')
    plt.plot(p, precision_list, label = 'Precision')
    plt.plot(p, recall_list, label = 'Recall')
    plt.plot(p, f1_score_list, label = 'F1 score')
    plt.plot(p, accuracy_list, label = 'Accuracy')
    plt.plot(p, mcc_list, label = 'MCC score')
    plt.legend(loc = 'lower left')
    plt.show()
    
if __name__ == '__main__':
    bcftools_vcf_file = "merged-normal.bam.mpileup.vcf.called.vcf"
    vcf_files = ["merged-normal.pileup.50.vcf", "merged-normal.pileup.55.vcf",
                 "merged-normal.pileup.60.vcf", "merged-normal.pileup.65.vcf",
                 "merged-normal.pileup.70.vcf", "merged-normal.pileup.75.vcf",
                 "merged-normal.pileup.80.vcf", "merged-normal.pileup.85.vcf",
                 "merged-normal.pileup.90.vcf", "merged-normal.pileup.95.vcf",
                 "merged-normal.pileup.100.vcf"]
    
    metrics(bcftools_vcf_file, vcf_files)
    
