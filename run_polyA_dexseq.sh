#!/bin/bash

supportTab="/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/ko_wt_support.tab"
code="threeprime_ko_wt"
outputDir="/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/dexseq/"${code}
inputDir="/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/"
species="mouse"

/share/apps/R-3.3.2/bin/Rscript polyA_dexseq.R --support.tab $supportTab --code $code --output.dir $outputDir --input.dir $inputDir --species $species
