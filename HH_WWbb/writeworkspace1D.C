#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TString.h"
using namespace RooFit ;

/* example to plot pdf in workspace
root -l file.root
workspace->cd()
RooPlot nnframe = NN->frame("title")
data_obs->plotOn(nnframe)
nnframe->Draw()
 */


// Signal: assume input signal cross section is 5pb from THistogram, rescale it into 1fb by 1e-3/5.0 before storing it into a workspace
// Inject 1fb signal into combine tool

static const unsigned int NProcess = 9;
static const unsigned int NSystematic = 9;
const TObjString processes[ NProcess ] = {"TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"};
const TObjString systematics[ NSystematic ] = {"CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale"};

void writeworkspace1D(int mass, char* channel, char* variable, float xmin, float xmax, char* infile, char* outfile){

    //TFile *file = new TFile("MTonly_HME_2D_xsec1pb_20180302_2D_NNbin10/Hhh_FinalBGYield_2Dlimits_2018-03-08_nnout_MTonly.root", "UPDATE");
    TFile *out = new TFile(outfile, "recreate");
    out->Close();
    
    char signalhist[50];
    sprintf(signalhist, "RadionM%d", mass);
    TObjString processes_histname[ NProcess ] = {"TT","sT","DY","data_untagged","TT_untagged","sT_untagged","ttV","VV", signalhist};
    TFile *file = new TFile(infile, "READ");
    if (!file->IsOpen()){
        cout <<"input file corrupted or not exist " << infile << endl;
        return;
    }

    //RooRealVar HME("HME","HME", 250.0, 1200.0);
    //RooRealVar NN(variable,variable, xmin, xmax);
    RooRealVar NN(variable, variable, xmin, xmax);
    cout <<"mass " << mass <<" variable "<< variable <<" xmin "<< xmin <<" xmax "<< xmax <<" infile "<< infile <<" outfile "<< outfile << endl;


    char wsname[50]; 
    sprintf(wsname, "%s_M%d", channel, mass);

    RooWorkspace *w = new RooWorkspace(wsname,"workspace") ;

    w->import(NN);
    //w->import(HME);


    TString channelname(channel);
    char data_obs_histname[50];
    sprintf(data_obs_histname, "data_obs_%s", wsname);
    cout << "data_obs_histname "<< data_obs_histname << endl;
    TH1F * hist_data = (TH1F*)  gDirectory->Get(data_obs_histname );
    RooDataHist hist_data_Roodata ("data_obs", "data_obs", RooArgList(NN), hist_data);
    //RooHistPdf hist_data_pdf ("data_obs", "data_obs", RooArgList(NN), hist_data_Roodata);

    w->import(hist_data_Roodata);
    //w->import(hist_data_pdf);
     
    for (unsigned int i=0; i < NProcess; i++){
    //for (unsigned int i=0; i < 1; i++){ //for test
        //skip DY for MuMu and ElEl channel
        cout <<" channel "<< channelname <<" process "<< processes[i].GetString() <<", start to import this process "<< endl;
        if ((channelname == "MuMu" || channelname == "ElEl") and processes[i].GetString().Contains("Drell_Yan")){
            cout <<" Skip !! "<< endl;
            continue;
        }
        if (channelname == "MuEl" and processes[i].GetString().Contains("untagged"))
            continue;
        
        //add nominal hist to workspace
        TString histname_nominal = processes_histname[i].GetString() +"_"+ channel; 
        TH1F * hist_nominal = (TH1F*)  gDirectory->Get( histname_nominal );
        if (processes[i].GetString().Contains("Signal"))
            hist_nominal->Scale(1e-3/5.0);
        //RooDataHist hist_nominal_Roodata (processes[i].GetString()+"_Roodata", processes[i].GetString()+"_Roodata", RooArgList(NN), hist_nominal);
        RooDataHist hist_nominal_Roodata (processes[i].GetString(), processes[i].GetString()+"_Roodata", RooArgList(NN), hist_nominal);
        //RooHistPdf hist_nominal_pdf (processes[i].GetString(), processes[i].GetString(), RooArgList(NN), hist_nominal_Roodata);

        w->import(hist_nominal_Roodata);
        //w->import(hist_nominal_pdf);
        if (processes[i].GetString().Contains("data"))
            continue;

        for (unsigned int j=0; j < NSystematic; j++){
        //for (unsigned int j=0; j < 1; j++){// for test
            TString histname_up = processes_histname[i].GetString() +"_"+ channel + "_"+ systematics[j].GetString()+"_up";
            TString histname_down = processes_histname[i].GetString() +"_"+ channel + "_"+ systematics[j].GetString()+"_down";
            TH1F * hist_up = (TH1F*)  gDirectory->Get( histname_up );
            TH1F * hist_down = (TH1F*)  gDirectory->Get( histname_down );
            //RooDataHist hist_up_Roodata (processes[i].GetString()+"_"+ systematics[j].GetString()+"_up_Roodata", processes[i].GetString()+"_"+ systematics[j].GetString()+"_up_Roodata", RooArgList(NN), hist_up);
            TString hist_sys_name = processes[i].GetString()+"_"+ systematics[j].GetString();

            //// QCD scale uncertainty for each process
            if (hist_sys_name.Contains("QCDscale") and hist_sys_name.Contains("untagged")){
                TString process_tmp = processes[i].GetString();
                hist_sys_name = processes[i].GetString()+"_"+ systematics[j].GetString() +process_tmp.ReplaceAll("_untagged","");
            }else if (hist_sys_name.Contains("QCDscale")  and !hist_sys_name.Contains("untagged"))
                hist_sys_name = processes[i].GetString()+"_"+ systematics[j].GetString() + processes[i].GetString();

            //// trigger eff uncertainty for each channel
            if (hist_sys_name.Contains("CMS_eff_trigger"))
                hist_sys_name = processes[i].GetString()+"_"+ systematics[j].GetString()+"_"+ channel;

            ///temperatory solution for mistake in signal PDF uncertainty
            if (processes[i].GetString().Contains("Signal") and systematics[j].GetString().Contains("CMS_pdf")){
                cout <<"process "<< hist_sys_name <<" up entries "<< hist_up->Integral() <<" down "<< hist_down->Integral();
            }

            RooDataHist hist_up_Roodata ( hist_sys_name + "Up", processes[i].GetString()+"_"+ systematics[j].GetString()+"_up_Roodata", RooArgList(NN), hist_up);
            //RooHistPdf hist_up_pdf (processes[i].GetString()+"_"+ systematics[j].GetString()+"Up", processes[i].GetString()+"_"+ systematics[j].GetString()+"_up", RooArgList(NN), hist_up_Roodata);
            //RooDataHist hist_down_Roodata (processes[i].GetString()+"_"+ systematics[j].GetString()+"_down_Roodata", processes[i].GetString()+"_"+ systematics[j].GetString()+"_down_Roodata", RooArgList(NN), hist_down);
            RooDataHist hist_down_Roodata (hist_sys_name +"Down", processes[i].GetString()+"_"+ systematics[j].GetString()+"_down_Roodata", RooArgList(NN), hist_down);
            //RooHistPdf hist_down_pdf (processes[i].GetString()+"_"+ systematics[j].GetString()+"Down", processes[i].GetString()+"_"+ systematics[j].GetString()+"_down", RooArgList(NN), hist_down_Roodata);
            w->import(hist_up_Roodata);
            //w->import(hist_up_pdf);
            w->import(hist_down_Roodata);
            //w->import(hist_down_pdf);
        }
    }

    w->Print();
    //w->writeToFile("rf_M400_workspace.root") ;
    w->writeToFile(outfile) ;
    //w->Write();


    file->Close();

    /*
    char histname_signal_MuMu[50];
    sprintf(histname_signal_MuMu, "signal_MuMu_M%d", mass);
    TH2F* hist_signal_MuMu = (TH2F*) gDirectory->Get(histname_signal_MuMu);
    //TH2F* hist_signal_MuMu = (TH2F*) gDirectory->Get("signal_MuMu_M400");
    RooDataHist signal_MuMu_data("signal_MuMu_data","signal_MuMu_data with (NN,HME)", RooArgList(NN,HME), hist_signal_MuMu);
    RooHistPdf signal_MuMu_pdf("signal_MuMu_pdf","signal_MuMu_pdf with (NN,HME)", RooArgList(NN,HME), signal_MuMu_data);

    char histname_TT_MuMu[50];
    sprintf(histname_TT_MuMu, "TT_MuMu_M%d", mass);
    TH2* hist_TT_MuMu = (TH2*) gDirectory->Get( histname_TT_MuMu );
    RooDataHist TT_MuMu_data("TT_MuMu_data","TT_MuMu_data with (NN,HME)", RooArgList(NN,HME), hist_TT_MuMu);
    RooHistPdf TT_MuMu_pdf("TT_MuMu_pdf","TT_MuMu_pdf with (NN,HME)", RooArgList(NN,HME), TT_MuMu_data);

    char histname_DY_MuMu[50];
    sprintf(histname_DY_MuMu, "DY_MuMu_M%d", mass);
    TH2* hist_DY_MuMu = (TH2*) gDirectory->Get(histname_DY_MuMu); 
    RooDataHist DY_MuMu_data("DY_MuMu_data","DY_MuMu_data with (NN,HME)", RooArgList(NN,HME), hist_DY_MuMu);
    RooHistPdf DY_MuMu_pdf("DY_MuMu_pdf","DY_MuMu_pdf with (NN,HME)", RooArgList(NN,HME), DY_MuMu_data);

    char histname_sT_MuMu[50];
    sprintf(histname_sT_MuMu, "sT_MuMu_M%d", mass);
    TH2* hist_sT_MuMu = (TH2*) gDirectory->Get(histname_sT_MuMu); 
    RooDataHist sT_MuMu_data("sT_MuMu_data","sT_MuMu_data with (NN,HME)", RooArgList(NN,HME), hist_sT_MuMu);
    RooHistPdf sT_MuMu_pdf("sT_MuMu_pdf","sT_MuMu_pdf with (NN,HME)", RooArgList(NN,HME), sT_MuMu_data);

    //Muel
    char histname_signal_MuEl[50];
    sprintf(histname_signal_MuEl, "signal_MuEl_M%d", mass);
    TH2*        hist_signal_MuEl = (TH2*) gDirectory->Get(histname_signal_MuEl);
    RooDataHist signal_MuEl_data("signal_MuEl_data","signal_MuEl_data with (NN,HME)", RooArgList(NN,HME), hist_signal_MuEl);
    RooHistPdf  signal_MuEl_pdf("signal_MuEl_pdf","signal_MuEl_pdf with (NN,HME)", RooArgList(NN,HME), signal_MuEl_data);

    char histname_TT_MuEl[50];
    sprintf(histname_TT_MuEl, "TT_MuEl_M%d", mass);
    TH2* hist_TT_MuEl = (TH2*) gDirectory->Get(histname_TT_MuEl);
    RooDataHist TT_MuEl_data("TT_MuEl_data","TT_MuEl_data with (NN,HME)", RooArgList(NN,HME), hist_TT_MuEl);
    RooHistPdf TT_MuEl_pdf("TT_MuEl_pdf","TT_MuEl_pdf with (NN,HME)", RooArgList(NN,HME), TT_MuEl_data);

    char histname_DY_MuEl[50];
    sprintf(histname_DY_MuEl, "DY_MuEl_M%d", mass);
    TH2* hist_DY_MuEl = (TH2*) gDirectory->Get(histname_DY_MuEl); 
    RooDataHist DY_MuEl_data("DY_MuEl_data","DY_MuEl_data with (NN,HME)", RooArgList(NN,HME), hist_DY_MuEl);
    RooHistPdf DY_MuEl_pdf("DY_MuEl_pdf","DY_MuEl_pdf with (NN,HME)", RooArgList(NN,HME), DY_MuEl_data);

    char histname_sT_MuEl[50];
    sprintf(histname_sT_MuEl, "sT_MuEl_M%d", mass);
    TH2* hist_sT_MuEl = (TH2*) gDirectory->Get(histname_sT_MuEl); 
    RooDataHist sT_MuEl_data("sT_MuEl_data","sT_MuEl_data with (NN,HME)", RooArgList(NN,HME), hist_sT_MuEl);
    RooHistPdf sT_MuEl_pdf("sT_MuEl_pdf","sT_MuEl_pdf with (NN,HME)", RooArgList(NN,HME), sT_MuEl_data);


    //Elel
    char histname_signal_ElEl[50];
    sprintf(histname_signal_ElEl, "signal_ElEl_M%d", mass);
    TH2*        hist_signal_ElEl = (TH2*) gDirectory->Get(histname_signal_ElEl);
    RooDataHist signal_ElEl_data("signal_ElEl_data","signal_ElEl_data with (NN,HME)", RooArgList(NN,HME), hist_signal_ElEl);
    RooHistPdf  signal_ElEl_pdf("signal_ElEl_pdf","signal_ElEl_pdf with (NN,HME)", RooArgList(NN,HME), signal_ElEl_data);

    char histname_TT_ElEl[50];
    sprintf(histname_TT_ElEl, "TT_ElEl_M%d", mass);
    TH2* hist_TT_ElEl = (TH2*) gDirectory->Get(histname_TT_ElEl);
    RooDataHist TT_ElEl_data("TT_ElEl_data","TT_ElEl_data with (NN,HME)", RooArgList(NN,HME), hist_TT_ElEl);
    RooHistPdf TT_ElEl_pdf("TT_ElEl_pdf","TT_ElEl_pdf with (NN,HME)", RooArgList(NN,HME), TT_ElEl_data);

    char histname_DY_ElEl[50];
    sprintf(histname_DY_ElEl, "DY_ElEl_M%d", mass);
    TH2* hist_DY_ElEl = (TH2*) gDirectory->Get(histname_DY_ElEl); 
    RooDataHist DY_ElEl_data("DY_ElEl_data","DY_ElEl_data with (NN,HME)", RooArgList(NN,HME), hist_DY_ElEl);
    RooHistPdf DY_ElEl_pdf("DY_ElEl_pdf","DY_ElEl_pdf with (NN,HME)", RooArgList(NN,HME), DY_ElEl_data);

    char histname_sT_ElEl[50];
    sprintf(histname_sT_ElEl, "sT_ElEl_M%d", mass);
    TH2* hist_sT_ElEl = (TH2*) gDirectory->Get(histname_sT_ElEl); 
    RooDataHist sT_ElEl_data("sT_ElEl_data","sT_ElEl_data with (NN,HME)", RooArgList(NN,HME), hist_sT_ElEl);
    RooHistPdf sT_ElEl_pdf("sT_ElEl_pdf","sT_ElEl_pdf with (NN,HME)", RooArgList(NN,HME), sT_ElEl_data);

    //RooWorkspace w("w",kTRUE);
    RooWorkspace *w = new RooWorkspace("w","workspace") ;

    w->import(NN);
    w->import(HME);

    w->import(signal_MuMu_data);
    w->import(TT_MuMu_data);
    w->import(DY_MuMu_data);
    w->import(sT_MuMu_data);

    w->import(signal_MuEl_data);
    w->import(TT_MuEl_data);
    w->import(DY_MuEl_data);
    w->import(sT_MuEl_data);

    w->import(signal_ElEl_data);
    w->import(TT_ElEl_data);
    w->import(DY_ElEl_data);
    w->import(sT_ElEl_data);

    w->import(signal_MuMu_pdf);
    w->import(TT_MuMu_pdf);
    w->import(DY_MuMu_pdf);
    w->import(sT_MuMu_pdf);

    w->import(signal_MuEl_pdf);
    w->import(TT_MuEl_pdf);
    w->import(DY_MuEl_pdf);
    w->import(sT_MuEl_pdf);

    w->import(signal_ElEl_pdf);
    w->import(TT_ElEl_pdf);
    w->import(DY_ElEl_pdf);
    w->import(sT_ElEl_pdf);
    w->Print();
    //w->writeToFile("rf_M400_workspace.root") ;
    w->writeToFile(outfile) ;
    //w->Write();


    file->Close();
    */

}
