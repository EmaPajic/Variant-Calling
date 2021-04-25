## Binomial Variant Caller from Pileup

This repo contains an implementation of a simple Variant Calling algorithm that assumes that nucleotide counts at one pileup position follow the same **Binomial** distribution, and uses **Maximum Likelihood Estimation (MLE)** for calling the variants.

[Youtube presentation](https://www.youtube.com/watch?v=TCuiXotgKk0&ab_channel=nikolaaleksic)

## Table of content

- [Binomial Variant Caller from Pileup](#binomial-variant-caller-from-pileup)
- [Table of content](#table-of-content)
- [How to use](#how-to-use)
- [What is Variant Calling?](#what-is-variant-calling)
- [Algorithm](#algorithm)
- [Results](#results)
- [Acknowlegments](#acknowlegments)
- [Licence](#licence)

## How to use

The algorithm takes an input file in the [Pileup format](https://en.wikipedia.org/wiki/Pileup_format), and produces a `.vcf` output file in the [Variant Call format](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

To run the tool:

* Install `pysam`, `numpy`, `scipy`. <span style="color:red">TODO: Create requirements.txt </span>
* Run `python main.py --input-file=PATH_TO_INPUT_FILE --output-file=PATH_TO_SAVE_OUTPUT --p=0.8`.

To see other possible parameters, call `python main.py --help`.

## What is Variant Calling?

Variant calling is the process of finding differences between a reference genome and an observed sample.
</br>
By comparing an observed genome with a well known reference genome, we can identify, store and query someone's genome easier.
</br>
Genomic variants can be simple, one example are single nucleotide variants (SNVs), and they can be complex, one example includes translocations.
<img src=images/genomic_variants.png>

Before identifying variants, we need to read the actual genome we want to analyze.
Unfortunately, there isn't a method to read the whole genome in one go, so DNA is read by fragmenting the genome into millions of molecules, and reading shorter [reads](https://en.wikipedia.org/wiki/Read_(biology)) instead. The reads are then aligned to reconstruct the whole genome.
</br>
</br>
The algorithm implemented in this repo identifies variants by looking at only one position at a time, for what the pileup format is perfect for.
</br>
[Pileup format](https://en.wikipedia.org/wiki/Pileup_format) gives us how many of each nucleotide **(A, G, C, T)** are at every position, and also how many insertion/deletition occurencies are at that position. It also gives us the information what the reference genome had, against what we will compare.

## Algorithm

After we have constrained the analysis to only one sequence position, we want to count how many of each **AGCT** letter and each **INDEL (insertion/deletition)** is at that position.
</br>
We will assume diploidy, hence there are at most two variants at one position. We will select between all letters and indels two with the largest count, remember their counts, and discard others.
</br>
</br>
To simplify, let's call them letters **a** and **b**, with counts **n<sub>a</sub>** and **n<sub>b</sub>**.
We will use the **Maximum Likelihood (MLE)** estimation to select the variants.
</br>
</br>
**P(variant|counts) = P(variant) * P(counts|variant) / P(counts)**
</br>
</br>
**MLE** ignores **P(variant)**, and chooses the most probable one by maximizing **P(counts|variant)**, *the likelihood function*.
</br>
</br>
All nucleotides at one position are either correct given the variant, or an error due to read extraction process. We will assume that all nucleotides are correct with probability **p**, and are independent of each other. There are three cases:

* Variant is **aa**, **n<sub>a</sub>** nucleotides are correct. Likelihood function is given by [Binomial distribution](https://en.wikipedia.org/wiki/Binomial_distribution) with parameters **n<sub>a</sub>** + **n<sub>b</sub>** and **p**, at position **n<sub>a</sub>**.
* Variant is **bb**, **n<sub>b</sub>** nucleotides are correct. Likelihood function is given by Binomial distribution with parameters **n<sub>a</sub>** + **n<sub>b</sub>** and **p**, at position **n<sub>b</sub>**.
* Variant is **ab**, all nucleotides are correct. We found that it is suboptimal to blend **a** and **b** into one correct letter, and have the likelihood function be equal to the same Binomial distribution at position **n<sub>a</sub>** + **n<sub>b</sub>**, because this option would be picked too often. To reduce that, we assumed that **a** and **b** should appear with equal probability **1/2**, and calculated the likelihood as probability that the counts are **n<sub>a</sub>** and **n<sub>b</sub>** in that setting.

After calculating the likelihoods, we will pick the most likely variant! That variant is compared to the referent genome, to determine the `genotype` and the `'alts'` field to store in the `VCF 4.2 format` output file. 
## Results

We used `merged-normal.bam` from [Seven Bridges Cancer Genomics Cloud ](https://www.cancergenomicscloud.org/) to test the algorithm. We produced `merged-normal.pileup` using `Bcftools Mpileup tool`, and compared our called variants with variants produces by `Bcftools Call tool`.
</br>
</br>
Having or not having a variant at one position in the `Bcftools Call tool` output will be our true/false labels, and comparing the output of our algorithm gives us `true positives`, `true negatives`, `false positives` and `false negatives`, which are used to calculate standard metrics.

We ran the algorithm for probabilities between 0.5 and 1, at every 0.05 increment, and the results are:
</br>
<img src=images/results.png>

The confusion matrix for **p** = 0.8 (best f1 score) is as follows:

<img src=images/confusion_matrix_0.8.png>

## Acknowlegments

This project was done for `Computational Genomics` course at *School of Electrical Engineering, University of Belgrade*.
</br>

You can check [the course github](https://github.com/vladimirkovacevic/gi-2021-etf) for useful learning materials.

`Bcftools Mpileup tool` and `Bcftools Call tool` were ran on [Seven Bridges Cancer Genomics Cloud ](https://www.cancergenomicscloud.org/), where we also downloaded the reference human genome and `merged-normal.bam` test file.


## Licence

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/EmaPajic/Variant-Calling/blob/main/LICENSE)
