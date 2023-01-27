#!/bin/bash

SCRAMARCH=slc7_amd64_gcc700
CMSSWVERSION=CMSSW_10_2_13
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAMARCH/cms/cmssw/$CMSSWVERSION/src ; eval `scramv1 runtime -sh` ; cd -

voms-proxy-init -voms cms
