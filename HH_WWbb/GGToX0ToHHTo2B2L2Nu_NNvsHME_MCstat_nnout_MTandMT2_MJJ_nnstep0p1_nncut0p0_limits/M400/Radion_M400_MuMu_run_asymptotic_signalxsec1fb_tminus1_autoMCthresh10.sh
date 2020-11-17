#!/bin/bash
echo 'start channel MuMu'
pushd /afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_nnout_MTandMT2_MJJ_nnstep0p1_nncut0p0_limits/M400/
# If workspace does not exist, create it once
if [ -f  GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu_combine_workspace.root ]; then
rm GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu_combine_workspace.root
fi

text2workspace.py GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu.dat -m 400 -o GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu_combine_workspace.root -P HiggsAnalysis.CombinedLimit.DYEstimationCombineModelTao:DYDataDrivenEstimationModelInstance
ValidateDatacards.py GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu.dat --mass 400 --printLevel 3 --jsonFile validateDatacard_M400_MuMu.json 
# Run limit

combine -M AsymptoticLimits -m 400 -n GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu_combine_workspace.root -s 1 -t -1 --verbose 1  &> GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuMu_signalxsec1fb_tminus1_autoMCthresh10.log
popd
echo 'finish channel MuMu'
