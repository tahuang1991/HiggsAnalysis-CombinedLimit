#! /bin/bash
pushd /afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_nnout_MTandMT2_MJJ_nnstep0p1_nncut0p0_limits/M400/
# If workspace does not exist, create it once
if [ ! -f  GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_combine_workspace.root ]; then
    exit 0 
fi

# Run impacts, -t -1 to fit on Asimov data

## example of likelihoodscan 
## $ root -l higgsCombineTest.MultiDimFit.mH125.root 
## limit->Draw("2*deltaNLL:r_ggH:r_qqH>>h(44,0,10,44,0,4)","2*deltaNLL<10","prof colz") 
combine --rMin 0 -m 400 -M MultiDimFit GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_combine_workspace.root -t -1 --verbose 3 &> GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_bestfit.log ## bestfit 
combine GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_combine_workspace.root -n M400_ElEl_rscan_expectS1 -M MultiDimFit --algo grid --points 2000 --rMin -1 --rMax 100 -m 400 --autoRange 1 --squareDistPoiStep -t -1 --fastScan  --expectSignal=1 

combineTool.py -M Impacts --rMax 100 -m 400 -n GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_S1 -d GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_combine_workspace.root --doInitialFit --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_analytic --expectSignal=1 -t -1 &> GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0_step1.log
combineTool.py -M Impacts --rMax 100 -m 400 -n GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_S1 -d GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_combine_workspace.root --doFits --parallel 4 --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_analytic --expectSignal=1 -t -1 &> GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0_step2.log
combineTool.py -M Impacts --rMax 100 -m 400 -n GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_S1 -d GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_combine_workspace.root -o GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0.json --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_analytic --expectSignal=1 -t -1 &> GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0_step3.log

plotImpacts.py -i GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0.json -o GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0 &> GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0_makeplot.log

cp GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M400_ElEl_impacts_signalxsec1fb_tminus1_peterExpectS0.pdf /afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_nnout_MTandMT2_MJJ_nnstep0p1_nncut0p0_limits/signalxsec1fb_tminus1_peterExpectS0_Impactplots/ 
popd
