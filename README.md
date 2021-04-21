## Binomial Variant Caller from Pileup

This repo contains an implementation of a simple Variant Calling algorithm that assumes that nucleotide counts at one pileup position follow the same **Binomial** distrubution, and uses **Maximum Likelihood Estimation (MLE)** for calling the variants.

[Youtube presentation](https://www.youtube.com/watch?v=TCuiXotgKk0&ab_channel=nikolaaleksic)

## Table of content

- [Binomial Variant Caller from Pileup](#binomial-variant-caller-from-pileup)
- [Table of content](#table-of-content)
- [How to use](#how-to-use)
- [What is Variant Calling?](#what-is-variant-calling)
- [Algorithm](#algorithm)
- [Results](#results)

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


## Algorithm

#TODO

## Results

#TODO
