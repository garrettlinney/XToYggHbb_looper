#!/bin/bash

cd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_30/ ; cmsenv ; cd - >/dev/null
source ~/miniconda3/etc/profile.d/conda.sh

INTERMEDIATEPKLFILES=$(python python/root_to_pickle.py $1)
conda activate fastparquet
python python/pickle_to_parquet.py $1
conda deactivate

rm $INTERMEDIATEPKLFILES
