#!/usr/bin/env bash

expName=$1
juicerCommand=$2
kindOfMatrix=$3
kindOfNormalization=$4
unitOfResolution=$5
resolution=$6
smoothing=$7
maxDistInteract=$8
regionsFileName=$9
hicFileName1="${10}"
hicFileName2="${11}"
outputFileName="${12}"

python hicCorr.py $expName $juicerCommand $kindOfMatrix $kindOfNormalization $unitOfResolution $resolution $smoothing $maxDistInteract $regionsFileName $hicFileName1 $hicFileName2 $outputFileName


