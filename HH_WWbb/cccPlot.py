import ROOT
import argparse
ROOT.gROOT.SetBatch(1)

poi = "r"
rMax = 10
filename="ccc.pdf"
infile="higgsCombineTest.ChannelCompatibilityCheck.mH125.root"
mass = 125
parser = argparse.ArgumentParser(description='runcccplot')
parser.add_argument("-p", "--POI", dest="poi", default="r", help="parameter of interesting . [Default: r]")
parser.add_argument("-rMax", "--rMax", dest="rMax", type=int, default=10, help="max R value . [Default: 10]")
parser.add_argument("-m", "--mH", dest="mH", type=int, default=125, help="mH  . [Default: 125]")
parser.add_argument("-o", "--output", dest="filename", default="ccc", help="output pdf name. [Default: ccc]")
parser.add_argument("-n", "--jobname", dest="jobname", default="ChannelCompatibilityCheck", help="job name, (as a part of inputfile). [Default: ChannelCompatibilityCheck]")
args = parser.parse_args()

poi = args.poi
rMax = args.rMax
mass = args.mH
filename = args.filename+".pdf"
infile = "higgsCombine"+args.jobname+".ChannelCompatibilityCheck.mH%d"%mass + ".root"


filein = ROOT.TFile(infile,"READ")
fit_nominal = filein.Get("fit_nominal");
fit_alternate = filein.Get("fit_alternate");
rFit = fit_nominal.floatParsFinal().find(poi);
print "posterfit rmin ",rFit.getMin()," rmax ",rFit.getMax()

prefix = "_ChannelCompatibilityCheck_"+poi+"_"

nChann = 0;
iter = fit_alternate.floatParsFinal().createIterator();

while True:
  a = iter.Next();
  if a==None: break
  if (prefix in a.GetName()): nChann+=1

frame = ROOT.TH2F("frame",";best fit #sigma/#sigma_{SM};",1,rFit.getMin(),min(rFit.getMax(),rMax),nChann,0,nChann)
    
iter.Reset(); 
iChann = 0; 
points = ROOT.TGraphAsymmErrors(nChann)

while True:
  a = iter.Next();
  if a==None: break
  if (prefix in a.GetName()):
    #RooRealVar *ri = (RooRealVar *) a;
    channel = a.GetName(); 
    channel.replace(prefix,"")
    points.SetPoint(iChann, a.getVal(), iChann+0.5);
    points.SetPointError(iChann, -a.getAsymErrorLo(), a.getAsymErrorHi(), 0, 0);
    print "iChann ",iChann," ",channel," value ",a.getVal()," error_low ",-a.getAsymErrorLo()," error high ",  a.getAsymErrorHi()
    iChann+=1;
    frame.GetYaxis().SetBinLabel(iChann, channel.replace(prefix,""));
    
c = ROOT.TCanvas("c","c",800,800)
c.cd()
c.SetTopMargin(0.07)
c.SetBottomMargin(0.12)
c.SetLeftMargin(0.12)

points.SetLineColor(ROOT.kRed);
points.SetLineWidth(3);
points.SetMarkerStyle(21);
frame.GetXaxis().SetTitleSize(0.05);
frame.GetXaxis().SetLabelSize(0.04);
frame.GetYaxis().SetLabelSize(0.06);
frame.Draw(); 
ROOT.gStyle.SetOptStat(0);
globalFitBand = ROOT.TBox(rFit.getVal()+rFit.getAsymErrorLo(), 0, rFit.getVal()+rFit.getAsymErrorHi(), nChann);
globalFitBand.SetFillColor(65);
globalFitBand.SetLineStyle(0);
globalFitBand.Draw("SAME");
globalFitLine = ROOT.TLine(rFit.getVal(), 0, rFit.getVal(), nChann);
globalFitLine.SetLineWidth(4);
globalFitLine.SetLineColor(214);
globalFitLine.Draw("SAME");
points.Draw("P SAME");

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.5*c.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right
latex2.DrawLatex(0.87, 0.95,"35.87 fb^{-1} (13 TeV)")
latex2.SetTextSize(0.7*c.GetTopMargin())
latex2.SetTextFont(62)
latex2.SetTextAlign(11) # align right
latex2.DrawLatex(0.15, 0.95, "CMS")
latex2.SetTextSize(0.6*c.GetTopMargin())
latex2.SetTextFont(52)
latex2.SetTextAlign(11)
latex2.DrawLatex(0.27, 0.95, "Tutorial")

ROOT.gPad.RedrawAxis()
c.SaveAs(filename);
