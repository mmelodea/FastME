//By Miquéias M. de Almeida 

#include "FastME.h"
#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#define cut 0.5					//Threshold to separate events (ideal cut 0.5 - MC #Sig and #Bkg equals)

using namespace std;

///****************************************************************************************************///
/// Calculate the distance between events using matrix with columms containing pT, eta, phi and lines  ///
/// containing (l1-l4,j1,j2) the leptons indice. The leptons Data and MC are organized by closure      ///
/// constraining that 2 leptons came from the same vertex (constraining Ressonance)      	       ///
///****************************************************************************************************///
TTree *FME::launch_FME(TString Final_State, TString out_name, TTree* Data_tree, TTree* MC_Sig_tree, TTree* MC_Bkg_tree, TString Branch_Name){
  int FS = -1;
  TString Model[3] = {"4l","4l2j","lv2j"};
       if(Final_State == Model[0]) FS = 0;
  else if(Final_State == Model[1]) FS = 1;
  else if(Final_State == Model[2]) FS = 2;
  else{
    cout<<"[Error] Model '"<<Final_State<<"' not defined!"<<endl;
    cout<<"Possible candidates are: "<<"4l, 4l2j ou lv2j"<<endl;
  }
    
  switch(FS){
    case 0:
      return FS4l_DR(out_name,Data_tree,MC_Sig_tree,MC_Bkg_tree,Branch_Name);
      break;
/*    
    case 1:
      return FS4l2j_DR(out_name,Data_tree,MC_Sig_tree,MC_Bkg_tree,Branch_Name);
      break;
      
    case 2:
      return FSlv2j_DR(TString out_name, TTree *Data_tree, TTree *MC_Sig_tree, TTree *MC_Bkg_tree, TString Branch_Name);
      break;
    
    default:
      return FS4l_DR(TString out_name, TTree *Data_tree, TTree *MC_Sig_tree, TTree *MC_Bkg_tree, TString Branch_Name);
      break;
*/  }
}

///Final State 4l only
TTree *FME::FS4l_DR(TString out_name, TTree* Data_tree, TTree* MC_Sig_tree, TTree* MC_Bkg_tree, TString Branch_Name){
  
  Float_t Data[4][3], MC[4][3];
  Data_tree->SetBranchAddress(Branch_Name,&Data);
  int ndata = Data_tree->GetEntries();
  MC_Sig_tree->SetBranchAddress(Branch_Name,&MC);
  int nsig = MC_Sig_tree->GetEntries();
  MC_Bkg_tree->SetBranchAddress(Branch_Name,&MC);
  int nbkg = MC_Bkg_tree->GetEntries();
  
  double minDR_toSig, minDR_toBkg, prob_sig_bkg, sig_frac, bkg_frac, min_dr_sig, min_dr_bkg;
  double pedestal = -99;
  int signal = 0, background = 0;
  TTree *FME_out = new TTree(out_name,"Fast Matrix Element Results");
  FME_out->SetDirectory(0);
  FME_out->Branch("minDR_toSig",&minDR_toSig);
  FME_out->Branch("minDR_toBkg",&minDR_toBkg);
  FME_out->Branch("P_SB",&prob_sig_bkg);
  FME_out->Branch("SigFrac",&sig_frac);
  FME_out->Branch("BkgFrac",&bkg_frac);
  
  for(int i=0; i<ndata; i++){
    if(i % (ndata/10) == 0) cout<<"Remaining Data: "<<ndata-i<<endl;
    Data_tree->GetEntry(i);
    
    //Reseting Variables
    minDR_toSig  = pedestal;
    minDR_toBkg  = pedestal;
    prob_sig_bkg = pedestal;
    sig_frac     = pedestal;
    bkg_frac	 = pedestal;
    min_dr_sig   = 1.E15;
    min_dr_bkg   = 1.E15;
    
    for(int s=0; s<nsig+nbkg; s++){
      if(s<nsig) MC_Sig_tree->GetEntry(s);
      else MC_Bkg_tree->GetEntry(s);
    
	double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
	double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
	int start = 0;
  
        do{
	    sum_dr1 = 0;
	    sum_dr2 = 0;
      
	    fdpt2    = pow(Data[start][0]-MC[start][0],2);//*Data[start][0][1];
	    fdeta2   = pow(Data[start][1]-MC[start][1],2);//*Data[start][1][1];
	    fdphi2   = pow(Data[start][2]-MC[start][2],2);//*Data[start][2][1];
	    sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
	    sdpt2    = pow(Data[start][0]-MC[start+1][0],2);//*Data[start][0][1];
	    sdeta2   = pow(Data[start][1]-MC[start+1][1],2);//*Data[start][1][1];
	    sdphi2   = pow(Data[start][2]-MC[start+1][2],2);//*Data[start][2][1];
	    sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
	  if(sum_dr1 < sum_dr2){
	    sum_dpt2  += fdpt2;
	    sum_deta2 += fdeta2;
	    sum_dphi2 += fdphi2;
	
	    sum_dpt2  += pow(Data[start+1][0]-MC[start+1][0],2);
	    sum_deta2 += pow(Data[start+1][1]-MC[start+1][1],2);
	    sum_dphi2 += pow(Data[start+1][2]-MC[start+1][2],2);
	  }
      
	  if(sum_dr2 < sum_dr1){
	    sum_dpt2  += sdpt2;
	    sum_deta2 += sdeta2;
	    sum_dphi2 += sdphi2;
	
	    sum_dpt2  += pow(Data[start+1][0]-MC[start][0],2);
	    sum_deta2 += pow(Data[start+1][1]-MC[start][1],2);
	    sum_dphi2 += pow(Data[start+1][2]-MC[start][2],2);
	  }

	  start += 2; 
	}while(start < 4);
	event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
	
	if(s<nsig)
	  if(event_distance < min_dr_sig) min_dr_sig = event_distance;
	  
	if(s>=nsig)
	  if(event_distance < min_dr_bkg) min_dr_bkg = event_distance;
    }    
    
      //Getting the discriminant value
      if(min_dr_sig == 0 && min_dr_bkg == 0){ 
	min_dr_sig = 1;
	min_dr_bkg = 1;
      }
      prob_sig_bkg = min_dr_sig/(min_dr_sig + min_dr_bkg);
      minDR_toSig = min_dr_sig;
      minDR_toBkg = min_dr_bkg;
      
      //Checking amount of kind choices
      if(prob_sig_bkg < cut) signal += 1;
      if(prob_sig_bkg > cut) background += 1;
      if(i == ndata-1){
	sig_frac = (signal/float(ndata))*100;
	bkg_frac = (background/float(ndata))*100;
	cout<<"sig = "<<sig_frac<<"%"<<endl;
	cout<<"bkg = "<<bkg_frac<<"%"<<endl;
      }
      
      //Saving Current Results
      FME_out->Fill();
  }
  //TFile *FME_Results = new TFile(name+".root","recreate");
  //FME_out->Write();
  //FME_Results->Close();
  return FME_out;
   
}

/*
///semi-leptonic case
TTree FME::FS_lv2j_DR(Float_t *Data, Float_t *MC){
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start=0;

    dpt  = (Data[start][0][0] - MC[start][0][0]);//*pt_res(start,MC[start][0][1],MC[start][1][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//*(1./MC[start][1][1]);
    dphi = (Data[start][2][0] - MC[start][2][0]);//*(1./MC[start][2][1]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    //MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0]);//*(1./MC[start][0][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//*(1./MC[start][1][1]);
    dphi = (Data[start][2][0] - MC[start][2][0]);//*(1./MC[start][2][1]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    //Jet 1
    start = 2;
    dpt  = (Data[start][0][0] - MC[start][0][0]);//*pt_res(start,MC[start][0][1],MC[start][1][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//MC[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0]);//MC[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    //Jet 2
    start = 3;
    dpt  = (Data[start][0][0] - MC[start][0][0]);//*pt_res(start,MC[start][0][1],MC[start][1][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//MC[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0]);//MC[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}*/
///====================================================================================================