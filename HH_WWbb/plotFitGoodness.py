import ROOT
import os
import sys

sys.argv.append( '-b' )


channellatex = {
    "MuMu":"#mu#mu",
    "MuEl":"e#mu",
    "ElEl":"ee",
    "MuMu_ElEl_MuEl":"combined"
}

def extranumber(s):
    num = []
    for t in s.split():
        try:
            num.append(float(t))
        except ValueError:
            pass

    return num


def extractTestStatisticfromtxtfile(logfile):
    logopen = open(logfile, "read")
    ts = []
    for line in logopen:
        if line.startswith("Best fit test statistic: "):
    	   ts = extranumber(line)
           break
    if len(ts) == 0:
      print("no best fit test statistic is found in logfile ", logfile)
      return 0
    else:
      #print("best fit test statistic is found in logfile ", ts[0])
      return ts[0]
            

def plotFitGoodness(mass, channel,plotdir, suffix):
    inputdir = "M%d"%mass
    treename = "limit"; branch = "limit"
    #higgsCombinefitgoodness_M400_MuMu.GoodnessOfFit.mH400.123456
    rootfile = os.path.join(inputdir, "higgsCombinefitgoodness_M{mass}_{ch}.GoodnessOfFit.mH{mass}.123456.root".format(mass=mass, ch = channel))
    tree = ROOT.TChain(treename)
    tree.Add(rootfile)
    #fitgoodness_M400_MuMu.log
    ts_data_f = os.path.join(inputdir, "fitgoodness_M{mass}_{ch}.log".format(mass=mass, ch = channel))
    ts_data = extractTestStatisticfromtxtfile(ts_data_f)
    c1 = ROOT.TCanvas("c1","c1",600, 800)
    hist = ROOT.TH1F("test_statistic","test_statistics",50, 10.0, 60.0)
    hist.SetTitle("Test statistics with Signal M=%d"%mass)
    hist.GetXaxis().SetTitle("Test_statistic")
    tree.Draw(branch+">> test_statistic")
    hmax = hist.GetMaximum()
    ts_data_line = ROOT.TLine(ts_data, 0.0, ts_data, hmax*0.7)
    ts_data_line.SetLineColor(ROOT.kRed)
    hist.Draw()
    ts_data_line.Draw("same")
    leg0 = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg0.SetFillColor(ROOT.kWhite)
    leg0.SetTextFont(42)
    #leg0.SetHeader()
    leg0.AddEntry(hist, "toys","l")
    leg0.AddEntry(ts_data_line, "data: %.1f"%ts_data,"l")
    leg0.Draw("same")

    tex1 = ROOT.TLatex(0.2,0.8, channellatex[channel])
    tex1.SetNDC(); tex1.SetTextSize(.045)
    tex1.Draw("same")


    c1.SaveAs(plotdir + "FitGoodness_M{mass}_{ch}_{suffix}.pdf".format(mass=mass, ch = channel, suffix=suffix))
    c1.SaveAs(plotdir + "FitGoodness_M{mass}_{ch}_{suffix}.C".format(mass=mass, ch = channel, suffix=suffix))

def plotAllFitGoodness(masslist, channels, plotdir, suffix):
    for mass in masslist:
      for ch in channels:
        plotFitGoodness(mass, ch, plotdir, suffix)



plotdir_fg = "FitGoodness_Plots/"
if not os.path.exists(plotdir_fg):
    os.system("mkdir "+plotdir_fg)
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
channels = ["ElEl","MuEl","MuMu", "MuMu_ElEl_MuEl"]   
#plotFitGoodness(400, "MuMu",plotdir_fg, "MTandMT2_MJJ_HME_antinncut0p8")
#plotFitGoodness(800, "MuMu",plotdir_fg, "MTandMT2_MJJ_HME_antinncut0p8")
plotAllFitGoodness(masslist, channels, plotdir_fg, "MTandMT2_MJJ_HME_antinncut0p8")



