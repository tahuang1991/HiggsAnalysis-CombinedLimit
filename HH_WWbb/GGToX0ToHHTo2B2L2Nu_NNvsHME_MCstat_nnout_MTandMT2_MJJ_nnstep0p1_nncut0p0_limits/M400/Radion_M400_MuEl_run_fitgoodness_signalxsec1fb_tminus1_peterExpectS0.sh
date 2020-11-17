#! /bin/bash
pushd /afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_nnout_MTandMT2_MJJ_nnstep0p1_nncut0p0_limits/M400/
combine -M GoodnessOfFit GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuEl.dat -m 400 --algo=saturated >& fitgoodness_M400_MuEl.log
combine -M GoodnessOfFit GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_MuEl.dat -m 400 --algo=saturated -t 1000  -n fitgoodness_M400_MuEl --expectSignal 0 
popd
