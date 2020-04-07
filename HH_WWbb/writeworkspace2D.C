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
#include "TH2D.h"
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

static const unsigned int NProcess = 9;
static const unsigned int NSystematic = 9;
const TObjString processes[ NProcess ] = {"TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"};
const TObjString systematics[ NSystematic ] = {"CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale"};

void writeworkspace2D(int mass, char* channel, char* variable, float xmin, float xmax, char* yvariable, float ymin, float ymax, char* infile, char* outfile){

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
    RooRealVar HME(yvariable, yvariable, ymin, ymax);
    cout <<"mass " << mass <<" variable "<< variable <<" xmin "<< xmin <<" xmax "<< xmax <<" yvariable "<< yvariable <<" ymin "<< ymin <<" ymax "<< ymax <<" infile "<< infile <<" outfile "<< outfile << endl;


    char wsname[50]; 
    sprintf(wsname, "%s_M%d", channel, mass);

    RooWorkspace *w = new RooWorkspace(wsname,"workspace") ;

    w->import(NN);
    w->import(HME);


    TString channelname(channel);
    char data_obs_histname[50];
    sprintf(data_obs_histname, "data_obs_%s", wsname);
    cout << "data_obs_histname "<< data_obs_histname << endl;
    TH2F * hist_data = (TH2F*)  gDirectory->Get(data_obs_histname );
    RooDataHist hist_data_Roodata ("data_obs", "data_obs", RooArgList(NN, HME), hist_data);
    //RooHistPdf hist_data_pdf ("data_obs", "data_obs", RooArgList(NN, HME), hist_data_Roodata);

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
        TH2F * hist_nominal = (TH2F*)  gDirectory->Get( histname_nominal );
        if (processes[i].GetString().Contains("Signal"))
            hist_nominal->Scale(1e-3/5.0);
        //RooDataHist hist_nominal_Roodata (processes[i].GetString()+"_Roodata", processes[i].GetString()+"_Roodata", RooArgList(NN), hist_nominal);
        RooDataHist hist_nominal_Roodata (processes[i].GetString(), processes[i].GetString()+"_Roodata", RooArgList(NN, HME), hist_nominal);
        //RooHistPdf hist_nominal_pdf (processes[i].GetString(), processes[i].GetString(), RooArgList(NN), hist_nominal_Roodata);

        w->import(hist_nominal_Roodata);
        //w->import(hist_nominal_pdf);
        if (processes[i].GetString().Contains("data"))
            continue;

        for (unsigned int j=0; j < NSystematic; j++){
        //for (unsigned int j=0; j < 1; j++){// for test
            TString histname_up = processes_histname[i].GetString() +"_"+ channel + "_"+ systematics[j].GetString()+"_up";
            TString histname_down = processes_histname[i].GetString() +"_"+ channel + "_"+ systematics[j].GetString()+"_down";
            TH2F * hist_up = (TH2F*)  gDirectory->Get( histname_up );
            TH2F * hist_down = (TH2F*)  gDirectory->Get( histname_down );
            //RooDataHist hist_up_Roodata (processes[i].GetString()+"_"+ systematics[j].GetString()+"_up_Roodata", processes[i].GetString()+"_"+ systematics[j].GetString()+"_up_Roodata", RooArgList(NN, HME), hist_up);
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

            RooDataHist hist_up_Roodata ( hist_sys_name + "Up", processes[i].GetString()+"_"+ systematics[j].GetString()+"_up_Roodata", RooArgList(NN, HME), hist_up);
            //RooHistPdf hist_up_pdf (processes[i].GetString()+"_"+ systematics[j].GetString()+"Up", processes[i].GetString()+"_"+ systematics[j].GetString()+"_up", RooArgList(NN, HME), hist_up_Roodata);
            //RooDataHist hist_down_Roodata (processes[i].GetString()+"_"+ systematics[j].GetString()+"_down_Roodata", processes[i].GetString()+"_"+ systematics[j].GetString()+"_down_Roodata", RooArgList(NN, HME), hist_down);
            RooDataHist hist_down_Roodata (hist_sys_name +"Down", processes[i].GetString()+"_"+ systematics[j].GetString()+"_down_Roodata", RooArgList(NN, HME), hist_down);
            //RooHistPdf hist_down_pdf (processes[i].GetString()+"_"+ systematics[j].GetString()+"Down", processes[i].GetString()+"_"+ systematics[j].GetString()+"_down", RooArgList(NN, HME), hist_down_Roodata);
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


}
