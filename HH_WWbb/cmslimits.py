import os
import ROOT
import re
import numpy as np
import datetime
import sys 
import csv
sys.argv.append( '-b' )
#or ROOT.gROOT.SetBatch(1)


import ROOT

def rescalesignal(rootfile, rdir, prefix, factor):
       rfile = ROOT.TFile(rootfile, "UPDATE")
       #rdir = "mjj_vs_NN_M400"
       a = rfile.Get(rdir)
       #print "a ",a, " ",a.GetListOfKeys()
       keys = a.GetListOfKeys()
       #print "keys ",keys
       for i, key in enumerate(keys):
	    #obj = ROOT.TH1F(key.ReadObj())
	    obj = key.ReadObj()
	    #print "obj ",obj
	    if prefix in obj.GetName():
		obj.Scale(factor)
		#print "hist name ", obj.GetName()
       rfile.Write("",ROOT.TObject.kOverwrite)

def rescalesignalall(masslist, workdir, factor):
    channels =   ["ElEl","MuEl","MuMu"]
    for mass in masslist:
       thisdir = workdir + "M%d.r7526/"%mass
       for ch in channels:
    	   rootfile = thisdir + "GGToX0ToHHTo2B2L2Nu_M%d_%s_shapes.root"%(mass, ch)
    	   newfile = thisdir + "GGToX0ToHHTo2B2L2Nu_M%d_%s_shapes_signalscale%d.root"%(mass, ch, int(factor))
    	   os.system("cp %s %s "%(rootfile, newfile))
    	   rdir = "mjj_vs_NN_M%d"%mass
    	   prefix = "ggX0HH%d"%mass
    	   rescalesignal(newfile, rdir, prefix, factor)
    

def extranumber(s):
    num = []
    for t in s.split():
        try:
            num.append(float(t))
        except ValueError:
            pass

    return num


def extractlimitfromtxtfile(logfile):
    logopen = open(logfile, "read")
    signal_xsec = 1.0
    #signal_xsec = 5000.0
    #percents = [2.5, 16.0, 50.0, 84.0, 97.5]
    limits = {}
    for line in logopen:
        #if line.startswith ("median expected limit: "):
    	if line.startswith("Expected 50.0%:"):
	    line = line.replace("Expected 50.0%:", "")
            nums = extranumber(line)
            limits[50.0]  = nums[0] * signal_xsec
        elif line.startswith("Expected  2.5%:"):
	    line = line.replace("Expected  2.5%:", "")
            nums = extranumber(line)
            limits[2.50] = nums[0] * signal_xsec
        elif line.startswith("Expected 16.0%:"):
	    line = line.replace("Expected 16.0%:", "")
            nums = extranumber(line)
            limits[16.0] = nums[0] * signal_xsec
        elif line.startswith("Expected 84.0%:"):
	    line = line.replace("Expected 84.0%:", "")
            nums = extranumber(line)
            limits[84.0] = nums[0] * signal_xsec
        elif line.startswith("Expected 97.5%:"):
	    line = line.replace("Expected 97.5%:", "")
	    nums = extranumber(line)
            limits[97.5] = nums[0] * signal_xsec
	elif line.startswith("Observed Limit:"):
	    nums = extranumber(line)
	    limits[-1] = nums[0] * signal_xsec
        else :
            pass
              
    print "limits ",limits
    return limits



def CombineLimitplots(filelist, histlist, masspoints, xtitle, legends, text, plotname):

    colors = [ROOT.kRed, ROOT.kMagenta+2,ROOT.kBlue+1, ROOT.kBlack, ROOT.kOrange]
    markers = [21,22,23, 20, 24]
    tfilelist = []
    for f in filelist:
        tfilelist.append( ROOT.TFile(f, "READ"))


    c1 = ROOT.TCanvas("c1","c1",600, 800)
    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    minx = min(masspoints)*0.8
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limits on production rate (fb)")
    b1.GetXaxis().SetTitle(xtitle)
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
    leg0.SetHeader( legends[-1] )
    #le0g.AddEntry(grdata,"observed","pl")

    drawUncertainty = True
    if drawUncertainty:
	i =  len(tfilelist)- 1
	tf = tfilelist[i]
        tf.cd()
	g_central = tf.Get(histlist[i]+"_central")
	g_central.SetLineColor(colors[i])
	g_central.SetMarkerColor(colors[i])
	g_central.SetMarkerStyle(markers[i])
	g_onesigma = tf.Get(histlist[i]+"_onesigma")
	g_twosigma = tf.Get(histlist[i]+"_twosigma")
	g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
	g_twosigma.SetLineStyle(2)
	g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
	g_onesigma.SetLineStyle(2)
	#leg0.AddEntry(g_data,"Observed","l")
	leg0.AddEntry(g_central,"Expected 95% upper limit","l")
	leg0.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
	leg0.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
	#g_central.Draw("lsame")
	g_twosigma.Draw("fe3same")
	g_onesigma.Draw("fe3same")


    leg = ROOT.TLegend(0.17,0.2,0.4,0.2+0.045*len(tfilelist))
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    leg.SetHeader("Observed")
    for i, tf in enumerate(tfilelist):
        tf.cd()
        g_data = tf.Get(histlist[i]+"_data")
        g_data.SetLineColor(colors[i])
        g_data.SetMarkerColor(colors[i])
        g_data.SetMarkerStyle(markers[i])
	g_central = tf.Get(histlist[i]+"_central")
	g_central.SetLineColor(colors[i])
	g_central.SetLineStyle(2)
	g_central.SetMarkerColor(colors[i])
	g_central.SetMarkerStyle(markers[i])
	g_central.Draw("lsame")
        g_data.Draw("lpsame")
        thisleg = leg.AddEntry(g_data, legends[i],"p")
        thisleg.SetTextColor(colors[i])
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    #tex1.Draw("same")
    
    leg0.Draw("same")
    leg.Draw("same")
    c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW_combined_v2.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW_combined_v2.C")

def makeBrazilPlot(masspoints, alllimits, xtitle, text, plotname):
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = []
    data = []
    for mass in masspoints:
        limits = alllimits[mass]	
        central.append(limits[50.0])
        twosigma_low.append(limits[2.5])
        onesigma_low.append(limits[16.0])
        onesigma_up.append(limits[84.0])
        twosigma_up.append(limits[97.5])
    	data.append(limits[-1])
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
    #print "allXpoints ",masspoints_all," onesigma ",onesigma_all," float masspoints ",masspoints_f
    g_data = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(data))
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
    g_data.SetLineWidth(2)
    g_data.SetLineStyle(1)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1)
    g_data.SetMarkerColor(ROOT.kBlack)

    #b1 = ROOT.TH1F("b2","b2",14, 250.0, 950.0)
    minx = min(masspoints)*0.9
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limits on production rate (fb)")
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
    g_data.Draw("lpsame")
    suffixname =  plotname.split("/")[-1]
    g_data.SetName("%s_data"%suffixname)
    g_central.SetName("%s_central"%suffixname)
    g_onesigma.SetName("%s_onesigma"%suffixname)
    g_twosigma.SetName("%s_twosigma"%suffixname)
    g_central.Write()
    g_data.Write()
    g_onesigma.Write()
    g_twosigma.Write()

    
    leg = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    leg.AddEntry(g_data,"Observed","pl")
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
    c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW.C")


def runImpacts(masslist, workdir, generatescripts, runscripts):
    channels =   ["ElEl","MuEl","MuMu", "MuMu_ElEl_MuEl"]
    scriptsuffix = "CMSHIG17006"
    #scriptsuffix = "ccc"
    channels = ["MuMu_ElEl_MuEl"]
    if generatescripts:
        #plotdir = workdir + "Impactsplots/"
        plotdir = workdir + "cccplots/"
	for mass in masslist:
	   thisdir = workdir + "M%d.r7526/"%mass
	   os.system("cp cccPlot.py "+thisdir)
	   for channel in channels:
	   	fname = thisdir + "Radion_M%d_%s_run_impacts_%s.sh"%(mass, channel, scriptsuffix)
		script = open(fname, "write")
		script.write("#! /bin/bash\n")
		script.write("pushd "+thisdir+"\n")
		script.write("# If workspace does not exist, create it once\n")
		script.write("if [ ! -f  GGToX0ToHHTo2B2L2Nu_M%d_%s_combine_workspace.root ]; then\n"%(mass, channel))
		script.write("    text2workspace.py GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}.dat -m {mass} -o GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root -P HiggsAnalysis.CombinedLimit.DYEstimationCombineModel:DYDataDrivenEstimationModelInstance\n".format(mass = mass, ch = channel))
		script.write("fi\n\n")
		script.write("# Run impacts\n\n")
		script.write("combineTool.py -M Impacts -m {mass} -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts -d GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root --doInitialFit &> GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}_step1.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))
		#script.write("combineTool.py -M Impacts -m {mass} -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts -d GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root --doFits --parallel 4 &> GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}_step2.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))
		script.write("combineTool.py -M Impacts -m {mass} -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts -d GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root --rMax 100 --doFits --parallel 4 &> /dev/null\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))##too much printout
		script.write("combineTool.py -M Impacts -m {mass} -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts -d GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root -o GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}.json &> GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}_step3.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))
		script.write("plotImpacts.py -i GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}.json -o GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix} &> GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}_makeplot.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))
		script.write("cp GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_impacts_{suffix}.pdf {plotdir} \n".format(mass = mass, ch = channel,  suffix =scriptsuffix, plotdir = plotdir))


		if channel == "MuMu_ElEl_MuEl":
		    ##copy ccPlot.py to local dir !!
		    script.write("combine -M ChannelCompatibilityCheck -m {mass} GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch} --saveFitResult --rMax 500\n".format(mass = mass, ch = channel))
		    script.write("python cccPlot.py  -m {mass} -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch} -o GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_ChannelCompatibilityCheck_{suffix}\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))
		    script.write("cp GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_ChannelCompatibilityCheck_{suffix}.pdf {plotdir} \n".format(mass = mass, ch = channel,  suffix =scriptsuffix, plotdir = plotdir))

		script.write("popd\n")

		os.system("chmod 775 "+fname)


    channels = ["MuMu_ElEl_MuEl"]
#SBATCH -o batchjobs_{mass}_{channel}-%A-%a.out
#SBATCH -e batchjobs_{mass}_{channel}-%A-%a.err
    if runscripts:
	for mass in masslist:
	   thisdir = workdir + "M%d.r7526/"%mass
	   for ich,channel in enumerate(channels):
	   	fname = thisdir + "Radion_M%d_%s_run_impacts_%s.sh"%(mass, channel, scriptsuffix)
		#os.system("source "+fname)
		#queue = "stakeholder"
		queue = "background-4g"
		jobscript = open("{0}/Send_Impacts_{1}_{2}_{3}.slrm".format(thisdir, mass, channel,  scriptsuffix), "w")
		jobscript.write("""#!/bin/bash
#SBATCH -J M{mass}{ich}
#SBATCH -p {partition}
#SBATCH -n1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=12:00:00
#SBATCH -o batchjobs_{mass}_{channel}_{suffix}.out
#SBATCH -e batchjobs_{mass}_{channel}_{suffix}.err
#SBATCH --ntasks-per-core=1

echo "starting at `date` on `hostname`"
echo "SLURM_JOBID=$SLURM_JOBID"
jobid=$SLURM_JOBID
source ~/.bashrc
source /cvmfs/cms.cern.ch/cmsset_default.sh
#cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/HH_WWbb/
eval `scramv1 runtime -sh`
source {filename}
echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format(mass = mass, ich = ich, channel = channel, filename = fname, partition=queue, suffix = scriptsuffix))
		jobscript.close()

		os.system("chmod 775 {0}/Send_Impacts_{1}_{2}_{3}.slrm".format(thisdir, mass, channel, scriptsuffix))
		os.system("sbatch {0}/Send_Impacts_{1}_{2}_{3}.slrm".format(thisdir, mass, channel, scriptsuffix))
    		#os.system("rm %s*step2*log"%thisdir)


def runCMSHIG17006Limits(masslist, workdir, generatescripts, runscripts):

    channels =   ["ElEl","MuEl","MuMu", "MuMu_ElEl_MuEl"]
    #channels =   ["ElEl","MuEl","MuMu"]
    chnames = ["ElEl","MuEl","MuMu", "all channels"]
    #scriptsuffix = "" ##HIG-17-006
    #scriptsuffix = "signalxsec10"
    scriptsuffix = "signalxsec1fb_fixedPDF"
    if generatescripts :
	for mass in masslist:
	    #thisdir = workdir + "M%d.r7526/"%mass
	    thisdir = workdir + "M%d/"%mass
	    for channel in ["ElEl","MuEl","MuMu", "MuMu_ElEl_MuEl"]:
		fname = thisdir+"Radion_M%d_%s_run_asymptotic_%s.sh"%(mass, channel, scriptsuffix)
		script = open(fname, "write")
		script.write("#! /bin/bash\n")
		script.write("pushd "+thisdir+"\n")
		script.write("# If workspace does not exist, create it once\n")
		#script.write("if [ ! -f  GGToX0ToHHTo2B2L2Nu_M%d_%s_combine_workspace.root ]; then\n"%(mass, channel))
		script.write("text2workspace.py GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}.dat -m {mass} -o GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root -P HiggsAnalysis.CombinedLimit.DYEstimationCombineModelTao:DYDataDrivenEstimationModelInstance\n".format(mass = mass, ch = channel))
		#script.write("fi\n\n")
		script.write("# Run limit\n\n")
		script.write("combine -M Asymptotic --rMax 500 -m {mass} -n GGToX0ToHHTo2B2L2Nu_M{mass}_{ch} GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_combine_workspace.root -S 1 &> GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_{suffix}.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix))
		script.write("popd\n")

		os.system("chmod 775 "+fname)

    if runscripts:
        alllimits = {}
	for channel in channels:
	    alllimits[channel] = {}
	for mass in masslist:
	    #thisdir = workdir + "M%d.r7526/"%mass
	    thisdir = workdir + "M%d/"%mass
	    limits = {}
	    for channel in channels:
		fname = thisdir+"Radion_M%d_%s_run_asymptotic_%s.sh"%(mass, channel, scriptsuffix)
	    	print "ch ",channel, " fname ",fname
		#os.system("source "+fname)
		#logfile = os.path.join(thisdir,  "GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}{suffix}.log".format(mass = mass, ch = channel, suffix =scriptsuffix)) ##HIG-17006
		logfile = os.path.join(thisdir,  "GGToX0ToHHTo2B2L2Nu_M{mass}_{ch}_{suffix}.log".format(mass = mass, ch = channel, suffix =scriptsuffix))
		limits = extractlimitfromtxtfile(logfile) 
    		print "limits ",limits
    		alllimits[channel][mass] = limits
	allplots = []
	rfiles = []
	histlist = []
	csvfile = open(workdir+'limits_%s.csv'%scriptsuffix, 'wb')
	writer = csv.writer(csvfile)
	writer.writerow([ ["channel_Mass"]  + [key for key in [-1, 2.5, 16.0, 50.0, 84.0, 97.5]] ])
	for i, channel in enumerate(channels):
	    for mass in masslist:
	    	rowname =  channel+"_%d"%mass
		#writer.writerow([ [rowname]  + [alllimits[channel][mass][key] for key in [-1, 2.5, 16.0, 50.0, 84.0, 97.5]] ])
	    	#for key, value in alllimits[channel][mass].items():
	    plotname = os.path.join(workdir, "Radion_"+channel+"_"+scriptsuffix)
	    makeBrazilPlot(masslist, alllimits[channel], "Radion Mass [GeV]", chnames[i], plotname)
	    allplots.append(plotname)
    	    rfiles.append(plotname+".root")
	    histlist.append("Radion_"+channel+"_"+scriptsuffix)
        rootfile =  os.path.join(workdir, "Radion_allchannels_"+scriptsuffix+".root")
	os.system("hadd -f "+rootfile+" "+' '.join(rfiles))
	allplotname = os.path.join(workdir, "Radion_final_"+scriptsuffix)
        #CombineLimitplots(rfiles, histlist, masslist, "Radion Mass [GeV]", chnames, "CMS-HIG-17-006", allplotname)
        CombineLimitplots(rfiles, histlist, masslist, "Radion Mass [GeV]", chnames, "HH#rightarrow bbWW #rightarrow bbl#nul#nu", allplotname)

	    


masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
#masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650]
#masslist = [750, 800, 900]

#masslist = [400]

#masslist = [260, 270]
#workdir = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/CMS_HIG_17_006/"
#workdir = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_limits/"
#workdir = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_limits_noPDF/"
#workdir = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_limits_20180704/"
workdir = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_2Dlimits_20180804/"
#workdir = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_2Dlimits_20180804_M750_800_900/"
#rescalesignalall(masslist, workdir, 1000.0)
runCMSHIG17006Limits(masslist, workdir, False, True)
#runImpacts(masslist, workdir, True, True)
