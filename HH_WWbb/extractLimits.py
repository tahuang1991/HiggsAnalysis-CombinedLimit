import os
import ROOT
import re
import numpy as np
import datetime
import Datacards
import sys


def writeintoworkspace_1D(masslist,inputdir, prefix_out, outdir):
    variable = "DNN"
    xmin = 0.12; xmax = 2.76
    for mass in masslist:
        outdir_mass = os.path.join(outdir, "M%d"%mass)
        infile = os.path.join(inputdir, "Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly_SignalM%d.root"%mass)
        os.system("mkdir -p "+outdir_mass)
        channels = ["MuMu","MuEl","ElEl"]
        for channel in channels:
            outfile = os.path.join(outdir_mass, prefix_out+"_M%d_%s_shapes.root"%(mass, channel))
            command = 'root -b -q "writeworkspace1D.C(%d, \\"%s\\", \\"%s\\", %f, %f,  \\"%s\\", \\"%s\\")"'%(mass, channel, variable, xmin, xmax, infile, outfile)
            print "writeintoworkspace command ",command
            os.system(command)
        ### example
        # root -b -q 'writeworkspace1D.C (400, "MuMu","NN", 0.0, 2.4, "HHbbWW_20180625_NNoutput_MjjCR_test/Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly_SignalM400.root","test.root")'
        ###

prefix_out = "GGToX0ToHHTo2B2L2Nu"

def generatedatacard(rootdir, masspoints, prefix, addobservation):
    #channels = ["MuMu","MuEl","ElEl"]
    channels = ["MuMu","ElEl","MuEl"]
    processes = ["TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"]
    datacard = Datacards.Datacards()
    for mass in masspoints:
        channels_rootfile = {}
        dir_mass = os.path.join(rootdir, "M%d"%mass)
        if not os.path.isdir(dir_mass):
            print "wrong directory for root files including shapes ", dir_mass
            sys.exit()
        for channel in channels:
            thisrootfile = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_shapes.root".format(mass = mass, channel = channel))
            channels_rootfile[channel] = thisrootfile
            if not os.path.isfile(thisrootfile):
                print "rootfile does not exist ", thisrootfile
                sys.exit()
            outfile = os.path.join(dir_mass, prefix+"_M{mass}_{channel}.dat".format(mass = mass, channel = channel))

            #print "channel ",channel
            datacard.generateDatacard( outfile, thisrootfile, channel, mass, addobservation)
        outfile_all = os.path.join(dir_mass, prefix+"_M{mass}_MuMu_ElEl_MuEl.dat".format(mass = mass))
        datacard.generateDatacard_multichannels(outfile_all, channels_rootfile, channels, mass, addobservation)

        



def alldatacards(rootdir, model, masspoints, addobservation, outdir, analysistype):
    #os.system("cp %s %s"%(rootfile, outdir))
    channels = ["MuMu","MuEl","ElEl"]
    processes = ["signal","TT","DY","sT"]
    for mass in masspoints:
        outfile = os.path.join(outdir, "%s_M%d.txt"%(analysistype, mass))
        suffix = "M%d"%mass
        i = 0
        #if mass<=400:
        #    i = 11
        #suffix = "RadionM%d_WP%i"%(mass, i) 
        rootfile = os.path.join(rootdir, "HME_%s_2D_bin10_M%d_workspace.root"%(model, mass))
        generatedatacard(outfile, rootfile, channels, processes, mass, suffix, addobservation, analysistype)



def alldatacards_HME(rootfile, masspoints, addobservation, outdir, analysistype):
    os.system("cp %s %s"%(rootfile, outdir))
    channels = ["MuMu","MuEl","ElEl"]
    processes = ["signal","TT","DY","sT"]
    #for mass in masspoints:
    for mass in masspoints:
      #for i in range(0, 12):
        i = 11
        if mass>400:
            i = 0
        suffix = "RadionM%d_WP%d"%(mass,i)
        outfile = os.path.join(outdir, "%s_%s.txt"%(analysistype, suffix))
        generatedatacard(outfile, rootfile, channels, processes, mass, suffix, addobservation, analysistype)


def extranumber(s):
    num = []
    for t in s.split():
        try:
            num.append(float(t))
        except ValueError:
            pass

    return num


def submitdatacards(cardname, method, mass):
    logfile = "expectedlimits_tmp.log"
    os.system("combine {cardname}  -M {method} -m {mass} -t -1 -s -1 > {logfile}".format(cardname = cardname, method = method, mass = mass, logfile = logfile))
    #os.system("combine {cardname}  -M {method} -t 100 -s -1 > {logfile}".format(cardname = cardname, method = method, logfile = logfile))
    logopen = open(logfile, "read")
    sigma_xsec = 1.0
    percents = [2.5, 16.0, 50.0, 84.0, 97.5]
    limits = {}
    for line in logopen:
        if not line.startswith("Expected"):
            continue
        #print "line ",line
        for percent in percents:
            pattern = "Expected%5.1f%%: r < "%(percent)
            #print " percent ",percent, " pattern ",pattern
            if re.match(pattern, line):
                limits[percent] = float(line.split(" ")[-1].rstrip()) * sigma_xsec
    #print "limits ",limits
    return limits



def submitdatacards_v2(cardname, method):
    logfile = "expectedlimits_tmp.log"
    os.system("combine {cardname}  -M {method} -t 100  -m {mass} -s -1 > {logfile}".format(cardname = cardname, method = method, mass = mass, logfile = logfile))
    logopen = open(logfile, "read")
    sigma_xsec = 1.0
    #percents = [2.5, 16.0, 50.0, 84.0, 97.5]
    limits = {}
    for line in logopen:
        if line.startswith ("median expected limit: "):
            nums = extranumber(line)
            limits[50.0]  = nums[0] * sigma_xsec
        elif line.startswith("   68% expected band :"):
            nums = extranumber(line)
            limits[16.0] = nums[0] * sigma_xsec
            limits[84.0] = nums[1] * sigma_xsec
        elif line.startswith("   95% expected band :"):
            nums = extranumber(line)
            limits[2.50] = nums[0] * sigma_xsec
            limits[97.50] = nums[1] * sigma_xsec
        else :
            pass
              
    print "limits ",limits
    return limits

def submitdatacards_v3(cardname, mass, method):
    logfile = "expectedlimits_tmp.log"
    outdir = "shapeM%d_%s"%(mass, method)
    os.system("mkdir -p "+outdir)
    os.system("combine {cardname}  -M {method} -m {mass} -t -1 -s -1 --plots --out {outdir} > {logfile}".format(cardname = cardname, method = method, mass = mass, outdir = outdir, logfile = logfile))

def runCombineTools(cardnameprefix, masspoints, method, logfile):
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        os.system("combine {cardname}  -M {method} -t 10000 -s -1 >> {logfile}".format(cardname = cardname, method = method, logfile = logfile))

def runImpacts(cardnameprefix, xpoints, logfile):
    for x in xpoints:
        cardname = cardnameprefix + "_M%d.txt"%x
        wsname = cardnameprefix + "_M%d.root"%x
        samplename = cardnameprefix + "_M%d"%x
        impacfile = samplename+"_impacts"
        os.system("text2workspace.py {cardname} -o {wsname} >> {logfile}".format(cardname = cardname, wsname = wsname, logfile = logfile))
        os.system("combineTool.py -M Impacts -d {wsname} -m 125 --doInitialFit >> {logfile}".format(wsname = wsname, logfile = logfile))
        os.system("combineTool.py -M Impacts -d {wsname} -m 125 --doFits --parallel 4 >> {logfile}".format(wsname = wsname, logfile = logfile))
        os.system("combineTool.py -M Impacts -d {wsname} -m 125 -o {impactfile}.json >> {logfile}".format(wsname = wsname, impactfile = impacfile,logfile = logfile))
        os.system("plotImpacts.py -i {impactfile}.json -o {impactfile} >> {logfile}".format(impactfile = impacfile, logfile = logfile))

def getlimits(logfile, mass):
    return True
    
def CombineLimitplots(filelist, histlist, masspoints, xtitle, legends, text, plotname):

    #colors = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta+2, ROOT.kOrange, ROOT.kBlack]
    #markers = [20, 20,  21,22,23,24]
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kBlack]
    markers = [20,  21,22,23,24]
    tfilelist = []
    for f in filelist:
        tfilelist.append( ROOT.TFile(f, "READ"))


    c1 = ROOT.TCanvas("c1","c1",700, 800)
    #c1 = ROOT.TCanvas("c1","c1",600, 800)

    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    minx = min(masspoints)*0.8
    maxx = max(masspoints)*1.2
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetXaxis().SetTitle(xtitle)
    b1.GetYaxis().SetTitle("95% CL limit on #sigma(gg#rightarrow X^{spin 0} #rightarrow HH) #times #it{Br}(HH #rightarrow b#bar{b}l#nul#nu) [fb]")
    b1.GetYaxis().SetTitleOffset(1.4)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #2000.0
    ylow = 1.0#1.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)
    b1.Draw()


    leg0 = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg0.SetFillColor(ROOT.kWhite)
    leg0.SetTextFont(42)
    #le0g.AddEntry(grdata,"observed","pl")

    #leg = ROOT.TLegend(0.52,0.48,0.9,0.52+0.05*len(tfilelist))
    leg = ROOT.TLegend(0.56,0.61,0.9,0.57+0.05*len(tfilelist))
    #leg.SetHeader("DNN inputs: Kinematic+")
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(.028)
    #leg.AddEntry(grdata,"observed","pl")
    for i, tf in enumerate(tfilelist):
        tf.cd()
        gcentralname = histlist[i]+"_central"
        if "data" in histlist[i]:
            gcentralname = histlist[i]
        g_central = tf.Get( gcentralname )
        g_central.SetLineColor(colors[i])
        g_central.SetMarkerColor(colors[i])
        g_central.SetMarkerStyle(markers[i])
        if i == 0:
            g_onesigma = tf.Get(histlist[i]+"_onesigma")
            g_twosigma = tf.Get(histlist[i]+"_twosigma")
            g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
            g_twosigma.SetLineStyle(2)
            g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
            g_onesigma.SetLineStyle(2)
            leg0.AddEntry(g_central,"Expected 95% upper limit","l")
            leg0.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
            leg0.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
            g_twosigma.Draw("fe3same")
            g_onesigma.Draw("fe3same")

        #if i ==0:
        #    g_central.Draw("lsame")
        #    thisleg = leg.AddEntry(g_central, legends[i],"l")
        #else:
        g_central.Draw("lpsame")
        thisleg = leg.AddEntry(g_central, legends[i],"pl")
        thisleg.SetTextColor(colors[i])
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    #tex1.Draw("same")
    
    leg0.Draw("same")
    leg.Draw("same")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_combined.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_combined.C")
    img = ROOT.TImage.Create();
    img.FromPad(c1)
    img.WriteImage(plotname+"_95Upperlmit_Nodata_combined.png")
    
def makeBrazilPlot(masspoints, central, onesigma_up, onesigma_low, twosigma_up, twosigma_low, xtitle, text, plotname):
    fakeerrors = [0.0]*len(masspoints)
    #g_onesigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(onesigma_low), np.array(onesigma_up))
    #g_twosigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(twosigma_low), np.array(twosigma_up))
    c1 = ROOT.TCanvas("c1","c1",600, 800)
    outfilename = plotname+".root"
    tfile = ROOT.TFile(outfilename,"UPDATE")
    
    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    onesigma_low.reverse()
    twosigma_low.reverse()
    onesigma_all = onesigma_up + onesigma_low
    twosigma_all = twosigma_up + twosigma_low
    masspoints_all = masspoints + list(reversed(masspoints))
    masspoints_f =  np.array(masspoints)+0.0
    print "allXpoints ",masspoints_all," onesigma ",onesigma_all," float masspoints ",masspoints_f
    g_central = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(central))
    g_onesigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(onesigma_all))
    g_twosigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(twosigma_all))


    #g_twosigma.SetFillColor(ROOT.kYellow)
    g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
    g_twosigma.SetLineStyle(2)
    g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
    g_onesigma.SetLineStyle(2)
    g_central.SetLineWidth(2)
    g_central.SetLineStyle(7)
    g_central.SetLineColor(9)
    g_central.SetMarkerStyle(20)
    g_central.SetMarkerSize(1)
    g_central.SetMarkerColor(9)

    #b1 = ROOT.TH1F("b2","b2",14, 250.0, 950.0)
    minx = min(masspoints)*0.8
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limit on #sigma(gg#rightarrow X #rightarrow HH) X #bf{#it{Br}}(HH #rightarrow b#bar{b}l#nul#nu) (fb)")
    b1.GetXaxis().SetTitle(xtitle)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #10000.0
    ylow = 1.0#10.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)

    b1.Draw()
    g_twosigma.Draw("fe3same")
    g_onesigma.Draw("fe3same")
    g_central.Draw("lpsame")
    samplename =  plotname.split("/")[-1]
    g_central.SetName("%s_central"%samplename)
    g_onesigma.SetName("%s_onesigma"%samplename)
    g_twosigma.SetName("%s_twosigma"%samplename)
    g_central.Write()
    g_onesigma.Write()
    g_twosigma.Write()

    
    leg = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    #leg.AddEntry(grdata,"observed","pl")
    leg.AddEntry(g_central,"Expected 95% upper limit","l")
    leg.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
    leg.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    tex1.Draw("same")
    
    leg.Draw("same")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_tmp.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_tmp.C")
    img = ROOT.TImage.Create();
    img.FromPad(c1)
    img.WriteImage(plotname+"_95Upperlmit_Nodata_tmp.png")
    #c1.SaveAs(plotname+"_95Upperlmit_Nodata_tmp.png")

def Limitplots(cardnameprefix, masspoints, method, text, plotname):
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = []
    method = "Asymptotic"
    """
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        submitdatacards_v3(cardname, mass, 'MaxLikelihoodFit')
    """ 
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        limits = submitdatacards(cardname, method, mass)
        print "cardname ",cardname," limits ",limits,"\n"
        central.append(limits[50.0])
        twosigma_low.append(limits[2.5])
        onesigma_low.append(limits[16.0])
        onesigma_up.append(limits[84.0])
        twosigma_up.append(limits[97.5])


def GetNLimits(cardnameprefix, mass, Ntimes, method, text, plotname):
    method = "Asymptotic"
    cardname = cardnameprefix + "_M%d.txt"%mass
    limits = submitdatacards(cardname, method, mass)
    hist = ROOT.TH1F("expectedlimits_M%d"%mass, "expectedlimits_M%d"%mass,100, limits[2.5], limits[97.5])
    for i in range(Ntimes):
        limits = submitdatacards_v2(cardname, method, mass)
        print "cardname ",cardname," limits ",limits,"\n"
        hist.Fill(limits[50.0])
        os.system("sleep 1")
    c1 = ROOT.TCanvas("c1","c1",800, 600)
    hist.Draw("hist")
    hist.SetLineWidth(2)
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    tex1.Draw("same")
    
    c1.SaveAs(plotname+"_95Upperlmit_centralexpected_Nodata_tmp.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_centralexpected_Nodata_tmp.C")
    


def Limitplots_HME(cardnameprefix, xpoints, mass,  method, text, plotname):
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = []
    method = "Asymptotic"
    """
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        submitdatacards_v3(cardname, mass, 'MaxLikelihoodFit')
    """ 
    for x in xpoints:
        suffix = "RadionM%d_WP%d"%(mass, x)
        cardname = cardnameprefix + "_%s.txt"%(suffix)
        limits = submitdatacards(cardname, method, mass)
        print "cardname ",cardname," limits ",limits,"\n"
        central.append(limits[50.0])
        twosigma_low.append(limits[2.5])
        onesigma_low.append(limits[16.0])
        onesigma_up.append(limits[84.0])
        twosigma_up.append(limits[97.5])

    xtitle = "WorkingPoint (NN cut)"
    makeBrazilPlot(xpoints, central, onesigma_up, onesigma_low, twosigma_up, twosigma_low, xtitle, text, plotname)



suffix = 'Radion_1D_nnoutMTandMJJ_xsec1pb'
output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
#datacarddir = "Datacards/"

#os.system("mkdir -p "+datacarddir)
inputdir = "HHbbWW_20180625_NNoutput_MjjCR_test/"
outdir = prefix_out+"_limits/"
os.system("mkdir -p "+outdir)
method = "Asymptotic"
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
#masslist = [400]
#writeintoworkspace_1D(masslist,inputdir, prefix_out, outdir)
#generatedatacard(outdir, masslist, prefix_out, True)
#alldatacards(rootfile, masslist, False, outdir)
cardnameprefix = os.path.join(outdir, "shape")
plotname = os.path.join(outdir, "Radion")
#masslist = [260, 500, 650,  900]
#Limitplots(cardnameprefix, masslist, method, plotname)
#modellist = ['MTonly','MT2only','MTandMT2','MTandMT2_MJJ','MTandMT2_HME','MTandMT2_HMEMJJ', 'MTandMJJ']
#modellist = ["MTonly","MTandMT2", "MTandMT2_MJJ"]
modellist = ["MTonly","MTandMT2"]
#modellist = ["MTonly"]
def getlimits_masspoint():
  gnames = []
  flist = []
  modelcombined = ""

  #analysistype = "cutandcount"
  analysistype = "shape_ws"

  for model in modellist:
    #suffix = 'Radion_2D_HME_NN_{model}_xsec1pb_{anatype}_tminus1_autoMCStats'.format(model = model, anatype = analysistype)
    #suffix = 'Radion_1D_HMEShape_{model}_xsec1pb_{anatype}_tminus1_autoMCStats'.format(model = model, anatype = analysistype)
    #output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
    #output_suffix = "2018-03-08_%s"%suffix
    #rootfile = "Hhh_FinalBGYield_xsec1pb_NN_nnout_{model}.root".format(model = model)
    #rootfile = "Hhh_FinalBGYield_xsec1pb_HMENNcut_nnout_{model}_2018-02-07.root".format(model = model)
    cwd = os.getcwd()
    datacarddir = os.path.join(cwd, "%s_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut035"%("MTandMT2"))
    #rootfile = os.path.join(datacarddir, "Hhh_FinalBGYield_2Dlimits_2018-03-13_nnout_%s.root"%model)
    rootdir = datacarddir
    print "rootfile ",rootdir
    
    #outdir = os.path.join(datacarddir, output_suffix)
    outdir = os.path.join(cwd, "%s_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut035_doublesyserr"%("MTandMT2"))

    #os.system("mkdir -p "+outdir)
    #if model != "MTonly":
    alldatacards(rootdir, model, masslist, False, outdir, analysistype)
    cardnameprefix = os.path.join(outdir, analysistype)
    print "cardnameprefix ",cardnameprefix

    plotname = os.path.join(outdir, "Radion_"+model+"_"+analysistype)
    #for mass in masslist:
    #for mass in [260, 500, 900]:
    #    GetNLimits(cardnameprefix, mass, 20, method, model+" M=%dGeV"%mass, plotname+"_M%d"%mass)
    
    #if model != "MTonly":
    Limitplots(cardnameprefix, masslist, method, model, plotname)
    flist.append(plotname + ".root")
    gnames.append("Radion_"+model+"_"+analysistype)
    modelcombined = modelcombined+model+"_"

  if len(modellist)>1:

      output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), modelcombined[:-1])
      outdir = os.path.join(datacarddir, output_suffix)
      os.system("mkdir -p "+outdir)
      #plotname = os.path.join(outdir, "Radion_"+modelcombined+analysistype)
      plotname = os.path.join(outdir, "Radion_"+modelcombined+"HME"+analysistype)
      legs = modellist
      if legs[-1] == "MTandMT2_MJJ":
          legs[-1] =  "MT,MT2andMJJ"
      CombineLimitplots(flist, gnames, masslist, "radion mass [GeV]", legs, " ", plotname)

#runImpacts("2018-02-07_Radion_1D_NN_MTonly_xsec1pb_AddObservation/shape", [270, 350, 450, 550, 600, 650, 750, 800], "runimpacts.log")
#getlimits_masspoint()
flist = []; gnames = [] ; legs = []
#flist = ["Radion_MuMu_ElEl_MuEl_CMSHIG17006.root","2018-02-07_Radion_1D_NN_MTonly_xsec1pb_tminus1/Radion_MTonly_NN.root","MTandMT2_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut035/Radion_MTonly_shape_ws.root", "MTandMT2_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut035_doublesyserr/Radion_MTonly_shape_ws.root"]
#gnames = ["Radion_MuMu_ElEl_MuEl","Radion_MTonly_NN", "Radion_MTonly_shape_ws","Radion_MTonly_shape_ws"]
#legs = ["HIG17006, Expected","DNN output shape","DNN output Vs HME(2D)","DNN output Vs HME(2D), double Sys."]
#plotname = "2018-03-09_MTonly_MTandMT2_MTandMT2_MJJ_2D/Radion_NNshape_HMEshape_NNVsHME2D_nodata_noHMEshape_doublesyserr"
flist= ["CMS_HIG_17_006/official/Radion_allchannels.root","GGToX0ToHHTo2B2L2Nu_limits_20180704/limitv1/Radion_allchannels_signalxsec1fb_fixedPDF.root", "GGToX0ToHHTo2B2L2Nu_2Dlimits_20180804/Radion_allchannels_signalxsec1fb_fixedPDF.root"]
legs = ["HIG17006","Tao:DNN output","Tao:DNN output vs HME"]
plotdir = "GGToX0ToHHTo2B2L2Nu_limits_2D_v2/"
#masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650]
for channel in ["MuMu","ElEl","MuEl"]:
#for channel in ["MuMu"]:
    plotname = os.path.join(plotdir, "2018-08-05_HIG17006_Tao_%s"%channel)
    gnames = ["Radion_%s_S1"%channel, "Radion_%s_signalxsec1fb_fixedPDF"%channel, "Radion_%s_signalxsec1fb_fixedPDF"%channel]
    CombineLimitplots(flist, gnames, masslist, "radion mass [GeV]", legs, channel, plotname)

#
#for nncut in ["03","035","04","045","05"]:
#    
#    cwd = os.getcwd()
#    datacarddir = os.path.join(cwd, "%s_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut%s"%("MTandMT2", nncut))
#    if nncut == "05":
#        datacarddir = os.path.join(cwd, "%s_HME_2D_xsec1fb_20180312_2D_NNbin10_nncut%s_v2"%("MTandMT2", nncut))
#    filename =  os.path.join(datacarddir, "Radion_MTonly_shape_ws.root")
#    flist.append( filename )
#    legs.append("DNN output >= 0.%s"%(nncut[1:]))
#    gnames.append("Radion_MTonly_shape_ws")
#
#plotname = "2018-03-13_MTonly_MTandMT2_MTandMT2_MJJ_2D_allNNcuts/Radion_NNshape_HMEshape_MTonly_NNcuts03_055"
#CombineLimitplots(flist, gnames, masslist, "radion mass [GeV]", legs, " ", plotname)
        
#CombineLimitplots(flist[1:], gnames[1:], masslist, "radion mass [GeV]", legs[1:], " ", plotname)


"""
def getlimits_WP():
  for model in modellist:
  #for mass in [260, 300, 350, 400, 500, 600, 750, 900]:
  #for mass in [300, 500, 900]:
   for mass in masslist:
    #suffix = 'Radion_1D_{model}_xsec1pb_M{mass}'.format(model = model, mass = mass)
    #Hhh_FinalBGYield_xsec1pb_NN_nnout_MTandMT2
    #rootfile = "Hhh_FinalBGYield_HME_nnout_{model}_HMETest_20180201.root".format(model = model)
    rootfile = "Hhh_FinalBGYield_xsec1pb_HMENNcut_nnout_{model}_2018-02-07.root".format(model = model)
    suffix = 'Radion_1D_HMEshape_{model}_xsec1pb_{anatype}_tminus1'.format(model = model, anatype = analysistype)
    output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
    outdir = os.path.join(datacarddir, output_suffix)
    os.system("mkdir -p "+outdir)
    alldatacards_HME(rootfile, [mass], False, outdir, "shape")
    cardnameprefix = os.path.join(outdir, "shape")

    text = model+" M=%d"%mass
    plotname = os.path.join(outdir, "Radion_"+model+"_M%d"%(mass))
    xpoints = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    Limitplots_HME(cardnameprefix, xpoints, mass,  method, text, plotname)


#getlimits_WP()

"""
