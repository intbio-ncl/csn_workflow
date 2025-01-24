# csn_workflow

This library provides the functions for carrying out the generation of co-evolutionary similarity networks. Coevolution matricies need to be created by CCMPred for a set of protein sequences, and protein sequences need to be aligned and then converted to psc format. 

It is recommended to use this fork of [CCMPred](https://github.com/mburridge96/CCMpred) for simple installation as master repo contains a compilation error. 

This is a refactored branch of the original csn_workflow that requires minimal setup and uses the latest libraries. 

## Requirements

To install this library, run:
```
pip install -r requirements.txt
```

To install CCMPred fork:
```
git clone --recursive https://github.com/mburridge96/CCMpred.git

cd CCMpred 

cmake . 

make
```

## Usage

To create Coevolution Similarity Networks, you first need to produce coevolution matrices by CCMPred for a set of protein sequences ([instructions here](https://github.com/soedinglab/CCMpred/wiki/FAQ)). You can download an already produced zipped folder of such matrices [here](https://data.ncl.ac.uk/articles/dataset/Trans241CoevMatrices_tar_gz/12555620). Decompress it.

Protein Multiple Sequence Alignments (MSAs) can be created using whichever software you desire, however these should be in clustal (.aln) format for network generation. These can also then be readily converted to .psc alignments which CCMPred requires as input using the `convert_alignment.py` script located in the `CCMPred/scripts` directory.

## Arguments

To run, you can either run each individual step in order:
    1. coev_net_creator
    2. create_aln_net
    3. computeCoev

Or, you can run the main.py file with the required arguments (file names will be defaulted)

## Example command

python3 main.py -n 360 -f output.mat -a proteins.aln

