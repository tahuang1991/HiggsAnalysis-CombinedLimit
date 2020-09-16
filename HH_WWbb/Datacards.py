import os
import ROOT
import re
import numpy as np
import datetime



class Datacards(object):
    """Class for generating datacards for higgs combine tool"""
    def __init__(self):
        """Initialize class """
        self.channels = ["MuMu", "MuEl","ElEl"]
        self.processes = ["Signal", "TTbar","SingleTop","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"]
        self.channels_processes = {
                "MuMu":["TTbar","SingleTop","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"],
                "ElEl":["TTbar","SingleTop","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"],
                "MuEl":["TTbar","SingleTop","Drell_Yan", "ttV","VV","Signal"],
                }
        ["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],
        self.SysDict_datacard = {
                "CMS_eff_b_heavy":{  "process":["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],### generally, should exclude "data_untagged"
                                     "description":"b-tagging uncertainty from b,c-jets", 
                                     "channel":self.channels,
                                     "isShape":True,
                                     },
                "CMS_eff_b_light":{  
                                     "process":["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],   
                                     "description":"b-tagging uncertainty from light-jets",
                                     "channel":self.channels,
                                     "isShape":True,
                                     },
                "CMS_pu":{ 
                           "process":["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],   
                           "description":"pielup uncertainty",
                           "channel":self.channels,
                           "isShape":True,
                           },
                "CMS_pdf":{
                            "process":["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],   
                            "description":"PDF uncertainty",
                            "channel":self.channels,
                            "isShape":True,
                            },
                "CMS_eff_trigger_MuMu":{
                                   "process":["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV"],   
                                   "description": "trigger eff uncertainty",
                                   "channel":["MuMu"],
                                   "isShape":True,
                                   },
                "CMS_eff_trigger_MuEl":{
                                   "process":["Signal", "TTbar","SingleTop","ttV","VV", "Drell_Yan"],   
                                   "description": "trigger eff uncertainty",
                                   "channel":["MuEl"],
                                   "isShape":True,
                                   },
                "CMS_eff_trigger_ElEl":{
                                   "process":["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV"],   
                                   "description": "trigger eff uncertainty",
                                   "channel":["ElEl"],
                                   "isShape":True,
                                   },
                "CMS_eff_e":{
                            "process": ["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],   
                            "description":"Electron eff",
                            "channel": ["MuEl","ElEl"],
                            "isShape":True,
                            },
                "CMS_eff_mu":{
                            "process": ["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],   
                            "description":"Muon eff",
                            "channel":["MuMu", "MuEl"],
                            "isShape":True,
                            },
                "CMS_iso_mu":{
                              "process": ["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],   
                              "description": "Muon Isolation",
                              "channel":["MuMu", "MuEl"],
                              "isShape":True,
                              },
                #"QCDscale":{
                #            "process": self.processes,
                #            "description":"QCD scale uncertainty",
                #            "channel":self.channels,
                #            "isShape":True,
                #            },
                "lumi_13TeV_2016":{
                            "process": ["Signal", "TTbar","SingleTop","TTbar_untagged","SingleTop_untagged", "ttV","VV", "Drell_Yan"],    
                            "description":"luminosity uncertainty",
                            "isShape":False,
                            "channel":self.channels,
                            "value": 1.025,
                        },
                "ttbar_xsec":{
                        "process": ["TTbar", "TTbar_untagged"],
                        "description":"TTbar cross section uncertainty",
                        "channel":self.channels,
                        "isShape":False,
                        "value": 1.053,
                        },
                "singleTop_xsec":{
                        "process":["SingleTop", "SingleTop_untagged"],
                        "description":"singletop cross section uncertainty",
                        "channel":self.channels,
                        "isShape":False,
                        "value":1.072,
                        },
                "dy_mc_xsec":{
                        "process":["Drell_Yan"],
                        "description":"Drell_Yan cross section uncertainty",
                        "channel":["MuEl"],
                        "isShape":False,
                        "value":1.05,
                        },
                "dy_rwgt_norm_MuMu":{
                        "process":["TTbar_untagged","SingleTop_untagged","data_untagged"],
                        "channel":["MuMu"],
                        "description":"Data-driven estimation normalization uncertainty",
                        "isShape":False,
                        "value":1.05,
                        },
                "dy_rwgt_norm_ElEl":{
                        "process":["TTbar_untagged","SingleTop_untagged","data_untagged"],
                        "channel":["ElEl"],
                        "description":"Data-driven estimation normalization uncertainty",
                        "isShape":False,
                        "value":1.05,
                        },
                }
        for process in ["TTbar","SingleTop","Drell_Yan", "ttV","VV","Signal"]:
            self.SysDict_datacard["QCDscale"+process] = {}
            self.SysDict_datacard["QCDscale"+process]["process"] = [process]
            self.SysDict_datacard["QCDscale"+process]["isShape"] = True
            self.SysDict_datacard["QCDscale"+process]["description"] = "QCD scale uncertainty for "+process
            if process in ["TTbar","SingleTop"]:
                self.SysDict_datacard["QCDscale"+process]["process"].append(process+"_untagged")
            if process != "Drell_Yan":
                self.SysDict_datacard["QCDscale"+process]["channel"] = self.channels
            else:
                self.SysDict_datacard["QCDscale"+process]["channel"] = ["MuEl"]

        
        #print "self.SysDict_datacard ",self.SysDict_datacard



    def readShapefromWorkspace(self, wsp, shapename, wsname, rootfile):
        hist = wsp.data(shapename)
        if not hist: hist = wsp.pdf(shapename)
        if not hist: hist = wsp.function(shapename)
        if not hist:
            raise RuntimeError, "Object %s in workspace %s in %s does not exist or it's neither a data nor a pdf" % (shapename, wsname, rootfile)
        return hist

    def getEntriesFromShape(self, wsname, shapename, rootfile, autoMC = False):
        rfile = ROOT.TFile(rootfile, "READ") ## read rate from histogram or workspace
        if autoMC:
          hist = rfile.Get(shapename)
          Entries = hist.Integral()
        if not autoMC:
          wsp = rfile.Get( wsname )
          hist = wsp.data(shapename)
          if not hist: hist = wsp.pdf(shapename)
          if not hist: hist = wsp.function(shapename)
          if not hist:
              raise RuntimeError, "Object %s in workspace %s in %s does not exist or it's neither a data nor a pdf" % (shapename, wsname, rootfile)
          Entries = hist.sumEntries()
          rfile.Close()
        return Entries


    def generateDatacard(self, outfile, rootfile, channel, mass, addobservation, autoMC = False):

        processes = self.channels_processes[channel]
        datacardfile = open(outfile, "write")
        jmax = len(processes)
        comment = "##datacard for GGToX0ToHHTo2B2L2Nu, Radion Model  Mass = %d GeV, %s"%(mass, channel)
        datacardfile.write(comment +"\n")
        datacardfile.write("imax 1\n")## only one here
        datacardfile.write("jmax %d\n"%(jmax-1))
        datacardfile.write("kmax *\n")  ## uncertainty sources
        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
        wsname = "{channel}_M{mass}".format(channel = channel, mass = mass)
        #datacardfile.write("shapes * * %s w:$PROCESS_$CHANNEL_pdf w:$PROCESS_$SYSTEMATIC_pdf\n"%(wsname))
        rootfile_local = rootfile.split("/")[-1]
        if autoMC: 
          datacardfile.write("shapes * * {rootfile} $PROCESS $PROCESS_$SYSTEMATIC\n".format(rootfile = rootfile_local))
        else: 
          datacardfile.write("shapes * * {rootfile} {wsname}:$PROCESS {wsname}:$PROCESS_$SYSTEMATIC\n".format(rootfile = rootfile_local, wsname = wsname))


        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
        datacardfile.write("bin\t\t"+wsname + "\n")
        #print "rootfile ",rootfile," used workspace to find hist integral"
        if addobservation:
            datacardfile.write("observation\t\t")
            oname = "data_obs"
            rate = self.getEntriesFromShape(wsname, oname, rootfile, autoMC)
	    datacardfile.write("{rate}\t".format(rate = rate))
            datacardfile.write("\n")
        else:
            datacardfile.write("#observation, not available yet\n")
        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
        nratelines = 4
        ratelineheads = ["bin","process","process","rate"]
        for iline,head in enumerate(ratelineheads):
            datacardfile.write("%-50s" % (head))
            for j, process in enumerate(self.channels_processes[channel]):
                if iline == 0:
                    datacardfile.write("%-20s" % (wsname))
                elif iline == 1:
                    datacardfile.write("%-20s" % (process))
                elif iline == 2:
                    if ("Signal" not in process) and "Radion" not in process: ## not signal
                        datacardfile.write("%-20d"%(j+1))
                    else:
                        datacardfile.write("%-20d"% (0))
                elif iline == 3:
                    oname = "{process}".format(process = process, channel = channel)
                    rate = self.getEntriesFromShape(wsname, oname, rootfile, autoMC)
                    datacardfile.write("%-20.3f"%(rate))
            datacardfile.write("\n")
        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out



        allSystematics = self.SysDict_datacard.keys()
        allSystematics.sort()
        for sys in allSystematics:
            if channel not in self.SysDict_datacard[sys]["channel"]:
                continue
            datacardfile.write("%-40s" % (sys))
            if self.SysDict_datacard[sys]["isShape"]:
                datacardfile.write("%-10s"%("shape"))
            else:
                datacardfile.write("%-10s"%("lnN"))
            if  channel not in self.SysDict_datacard[sys]["channel"]:
                continue
            for process in self.channels_processes[channel]:
                if channel in self.SysDict_datacard[sys]["channel"] and process in self.SysDict_datacard[sys]["process"]:
                    if self.SysDict_datacard[sys]["isShape"]:
                        datacardfile.write("%-20.1f"%(1.0))
                    else:
                        datacardfile.write("%-20.3f"%self.SysDict_datacard[sys]["value"])
                else:
                    datacardfile.write("%-20s"% ("-"))
            datacardfile.write("\n")

   
        if autoMC: 
          #datacardfile.write("* autoMCStats 0.0 0 1\n") 
          datacardfile.write("* autoMCStats 10 0 1 \n") 
        datacardfile.close()


    def generateDatacard_multichannels(self, outfile, channels_rootfile, channels, mass, addobservation, autoMC = False):


        datacardfile = open(outfile, "write")
        #jmax = len(self.channels_processes["MuMu"]) + 1
        jmax = len(self.processes)
        comment = "##datacard for GGToX0ToHHTo2B2L2Nu, Radion Model  Mass = %d GeV, "%(mass) + '-'.join(channels)
        datacardfile.write(comment +"\n")
        datacardfile.write("imax %d\n"%(len(channels)))## only one here
        datacardfile.write("jmax %d\n"%(jmax-1))
        datacardfile.write("kmax *\n")  ## uncertainty sources
        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
        for channel in channels:
            rootfile = channels_rootfile[channel]
            rootfile_local = rootfile.split("/")[-1]
            wsname = "{channel}_M{mass}".format(channel = channel, mass = mass)
            #datacardfile.write("shapes * * %s w:$PROCESS_$CHANNEL_pdf w:$PROCESS_$SYSTEMATIC_pdf\n"%(wsname))
            if autoMC: 
              datacardfile.write("shapes * {channel} {rootfile} $PROCESS $PROCESS_$SYSTEMATIC\n".format(rootfile = rootfile_local, wsname = wsname, channel = channel))
            else: 
              datacardfile.write("shapes * {channel} {rootfile} {wsname}:$PROCESS {wsname}:$PROCESS_$SYSTEMATIC\n".format(rootfile = rootfile_local, wsname = wsname, channel = channel))

        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
        
        datacardfile.write("bin\t")
        for channel in channels:
            datacardfile.write("\t"+ channel + "\t")
            #print "rootfile ",rootfile," used workspace to find hist integral"
        datacardfile.write("\n")

        if addobservation:
            datacardfile.write("observation ")
            for channel in channels:
                wsname = "{channel}_M{mass}".format(channel = channel, mass = mass)
                rootfile = channels_rootfile[channel]
                oname = "data_obs"
                #hist = self.readShapefromWorkspace( wsp, oname, wsname, rootfile)
                #rate = hist.sumEntries()
                rate = self.getEntriesFromShape(wsname, oname, rootfile, autoMC)
                datacardfile.write("\t{rate}\t".format(rate = rate))
            datacardfile.write("\n")
        else:
            datacardfile.write("#observation, not available yet\n")

        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out
        nratelines = 4
        ratelineheads = ["bin","process","process","rate"]

        for iline,head in enumerate(ratelineheads):
            datacardfile.write("%-50s" % (head))
            for channel in channels:
                rootfile = channels_rootfile[channel]
                wsname = "{channel}_M{mass}".format(channel = channel, mass = mass)
                for j, process in enumerate(self.channels_processes[channel]): 
                    if iline == 0:
                        datacardfile.write("%-20s" % (channel))
                    elif iline == 1:
                        datacardfile.write("%-20s" % (process))
                    elif iline == 2:
                        datacardfile.write("%-20d"%(self.processes.index(process)) )
                    elif iline == 3:
                        oname = "{process}".format(process = process, channel = channel)
                        rate = self.getEntriesFromShape(wsname, oname, rootfile, autoMC)
                        datacardfile.write("%-20.3f"%(rate))
            datacardfile.write("\n")
        datacardfile.write("#-----------------------------------------------------------------------------\n")## comment out



        allSystematics = self.SysDict_datacard.keys()
        allSystematics.sort()
        for sys in allSystematics:
            datacardfile.write("%-40s" % (sys))
            if self.SysDict_datacard[sys]["isShape"]:
                datacardfile.write("%-10s"%("shape"))
            else:
                datacardfile.write("%-10s"%("lnN"))
            for channel in channels:
                for process in self.channels_processes[channel]:
                    if channel in self.SysDict_datacard[sys]["channel"] and process in self.SysDict_datacard[sys]["process"]:
                        if self.SysDict_datacard[sys]["isShape"]:
                            datacardfile.write("%-20.1f"%(1.0))
                        else:
                            datacardfile.write("%-20.3f"%self.SysDict_datacard[sys]["value"])
                    else:
                        datacardfile.write("%-20s"% ("-"))
            datacardfile.write("\n")


        if autoMC: 
          #datacardfile.write("* autoMCStats 0.0 0 1\n")    
          datacardfile.write("* autoMCStats 10 0 1 \n")    
        datacardfile.close()




