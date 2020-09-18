import ROOT
import sys
import os

sys.argv.append( '-b' )


channellatex = {
    "MuMu":"#mu#mu",
    "MuEl":"e#mu",
    "ElEl":"ee",
    "MuMu_ElEl_MuEl":"combined"
}

def plotLikelihoodScan1D(mass, channel, variable, xmin, xmax, plotdir, suffix):
    treename = "limit"

    infile = "higgsCombineGGToX0ToHHTo2B2L2Nu_NNvsHME_MCstat_M800_MuMu_likelihoodscan_r_expectS1.MultiDimFit.mH800.root"

    tree = ROOT.TChain(treename)
    tree.Add(infile)

    hist = ROOT.TProfile("hist", "likelihood scan", 100, xmin, xmax)

    #tree.Draw("2*deltaNLL:"+variable+" >> hist", "2*deltaNLL<10")
    tree.Draw("2*deltaNLL:"+variable+" >> hist")
    hist.Print("ALL")


    c1 = ROOT.TCanvas("c1","c1",600, 800)

    c1.SetGridx()
    c1.SetGridy()                                                                                                                      
    c1.SetTickx()
    c1.SetTicky()

    hist.Draw()
    hist.GetXaxis().SetTitle(variable)
    hist.GetYaxis().SetTitle("2*deltaNLL")
    hist.SetTitle("Likelihood scan with Signal M=%d"%mass)

    tex1 = ROOT.TLatex(0.2,0.8, channellatex[channel])
    tex1.SetNDC(); tex1.SetTextSize(.045)
    tex1.Draw("same")

    c1.SaveAs(plotdir+"Likelihoodscan_2deltaNLL_{var}_M{mass}_{ch}_{suffix}.pdf".format(mass=mass, ch = channel, suffix=suffix,var=variable))


plotLikelihoodScan1D(800, "MuMu", "r", -1.0, 8.0, "./", "MTandMT2_MJJ_NNvsHME")
