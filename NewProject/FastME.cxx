///FastME Constructor
///Author: Miquéias M. de Almeida 

#include "FastME.h"
#include <iostream>
#include <ctime>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#define pedestal -99				///Reset Value to Variables
#define cut 0.5					///Threshold to Separate Events (Ideal Cut 0.5 - MC #Sig and #Bkg Equals)

using namespace std;

///:::::::: PHASE ESPACE DATA-MC DISTANCE COMPUTER FOR EQUAL RESSONANCE COMPONENTS :::::::::
double std_DR(const int comp, Float_t Data[comp][3][2], Float_t MC[comp][3][2]){
  double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  do{
      sum_dr1 = 0;
      sum_dr2 = 0;
      
      fdpt2    = pow(Data[start][0][0]-MC[start][0][0],2)/Data[start][0][1];
      fdeta2   = pow(Data[start][1][0]-MC[start][1][0],2)/Data[start][1][1];
      fdphi2   = pow(Data[start][2][0]-MC[start][2][0],2)/Data[start][2][1];
      sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
      sdpt2    = pow(Data[start][0][0]-MC[start+1][0][0],2)/Data[start][0][1];
      sdeta2   = pow(Data[start][1][0]-MC[start+1][1][0],2)/Data[start][1][1];
      sdphi2   = pow(Data[start][2][0]-MC[start+1][2][0],2)/Data[start][2][1];
      sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
    if(sum_dr1 < sum_dr2){
      sum_dpt2  += fdpt2;
      sum_deta2 += fdeta2;
      sum_dphi2 += fdphi2;
  
      sum_dpt2  += pow(Data[start+1][0][0]-MC[start+1][0][0],2)/Data[start+1][0][1];
      sum_deta2 += pow(Data[start+1][1][0]-MC[start+1][1][0],2)/Data[start+1][1][1];
      sum_dphi2 += pow(Data[start+1][2][0]-MC[start+1][2][0],2)/Data[start+1][2][1];
    }
      
    if(sum_dr2 < sum_dr1){
      sum_dpt2  += sdpt2;
      sum_deta2 += sdeta2;
      sum_dphi2 += sdphi2;

      sum_dpt2  += pow(Data[start+1][0][0]-MC[start][0][0],2)/Data[start+1][0][1];
      sum_deta2 += pow(Data[start+1][1][0]-MC[start][1][0],2)/Data[start+1][1][1];
      sum_dphi2 += pow(Data[start+1][2][0]-MC[start][2][0],2)/Data[start+1][2][1];
    }

    start += 2; 
  }while(start < comp);
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  
  return event_distance;
}
///==========================================================================================







///::::::::::::::::::::::::::::::::::	FAST MATRIX ELEMENT CLASS MEMBERS DEFINITIONS	::::::::::::::::::::::::::::::::::
int FME::launchFME(TString Final_State, TString Out_Name, TString Data_Path, TString MC_Sig_Path, TString MC_Bkg_Path, TString Tree_Name, TString Branch_Name){
  const int nModel = 3;
  int FS = -1;
  TString Model[nModel] = {"4l","4l2j","lv2j"};
       if(Final_State == Model[0]) FS = 0;
  else if(Final_State == Model[1]) FS = 1;
  else if(Final_State == Model[2]) FS = 2;
  else{
    cout<<"[Error] Model '"<<Final_State<<"' not defined!"<<endl;
    cout<<"Possible candidates are: ";
    for(int nM=0; nM<nModel; nM++) cout<<Model[nM]<<", "<<endl;
  }
  
  ///Preparing Inputs
  cout<<"Preparing Inputs..."<<endl;
  TFile *fData = TFile::Open(Data_Path);
  TFile *fMC_Sig = TFile::Open(MC_Sig_Path);
  TFile *fMC_Bkg = TFile::Open(MC_Bkg_Path);
  TTree *Data_Tree = (TTree*)fData->Get(Tree_Name);
  TTree *MC_Sig_Tree = (TTree*)fMC_Sig->Get(Tree_Name);
  TTree *MC_Bkg_Tree = (TTree*)fMC_Bkg->Get(Tree_Name);
    
  switch(FS){
    case 0:
      return FS_4l(Out_Name,Data_Tree,MC_Sig_Tree,MC_Bkg_Tree,Branch_Name);
      break;
    
    case 1:
      return FS_4l2j(Out_Name,Data_Tree,MC_Sig_Tree,MC_Bkg_Tree,Branch_Name);
      break;
      
    case 2:
      return FS_lv2j_DR(TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
      break;
    
    default:
      return FS_4l(TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
      break;
  }
}


///:::::::::::::::::::::	FULL LEPTONIC CASE	   :::::::::::::::::::::
int FME::FS_4l(TString Out_Name, TTree* Data_Tree, TTree* MC_Sig_Tree, TTree* MC_Bkg_Tree, TString Branch_Name){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
  
  cout<<"\n:::::::::::::::: Starting FastME Processing ::::::::::::::::"<<endl;
  cout<<"Model Chosen:  4l"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  Float_t Data[4][3][2], MC[4][3][2];
  Data_Tree->SetBranchAddress(Branch_Name,&Data);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Branch_Name,&MC);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Branch_Name,&MC);
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  double minDR_toSig, minDR_toBkg, prob_sig_bkg, sig_frac, event_weight;
  double bkg_frac, min_dr_sig, min_dr_bkg, sig_event_weight, bkg_event_weight;
  int signal = 0, background = 0, sig_neighbors, bkg_neighbors;
  TTree *FME_out = new TTree(Out_Name,"Fast Matrix Element Results");
  FME_out->SetDirectory(0);
  FME_out->Branch("minDR_toSig",&minDR_toSig);
  FME_out->Branch("minDR_toBkg",&minDR_toBkg);
  FME_out->Branch("SigFrac",&sig_frac);
  FME_out->Branch("BkgFrac",&bkg_frac);
  FME_out->Branch("P_SB",&prob_sig_bkg);
  FME_out->Branch("WSig_ToEvent",&sig_event_weight);
  FME_out->Branch("WBkg_ToEvent",&bkg_event_weight);
  FME_out->Branch("Event_Weight",&event_weight);

  
  for(int i=0; i<ndata; i++){
    if(i % (ndata/10) == 0) cout<<"Remaining Data: "<<ndata-i<<endl;
    Data_Tree->GetEntry(i);
    
    //Reseting Variables
    minDR_toSig   	= pedestal;
    minDR_toBkg   	= pedestal;
    prob_sig_bkg  	= pedestal;
    sig_frac      	= pedestal;
    bkg_frac	  	= pedestal;
    sig_event_weight 	= pedestal;
    bkg_event_weight 	= pedestal;
    min_dr_sig    	= 1.E15;
    min_dr_bkg    	= 1.E15;
    sig_neighbors 	= 0;
    bkg_neighbors 	= 0;
    
    for(int s=0; s<nsig; s++){
      MC_Sig_Tree->GetEntry(s);
      double dr_test = std_DR(4,Data,MC);
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
      if(dr_test < phs_radius) sig_neighbors += 1;
    }
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      double dr_test = std_DR(4,Data,MC);
      if(dr_test < min_dr_bkg) min_dr_bkg = dr_test;
      if(dr_test < phs_radius) bkg_neighbors += 1;
    }    
    
    ///Getting the discriminant value
    if(min_dr_sig == 0 && min_dr_bkg == 0){ 
      min_dr_sig = 1;
      min_dr_bkg = 1;
    }
      prob_sig_bkg = min_dr_sig/(min_dr_sig + min_dr_bkg);
      minDR_toSig = min_dr_sig;
      minDR_toBkg = min_dr_bkg;
      
    ///Checking amount of kind choices
    if(prob_sig_bkg < cut) signal += 1;
    if(prob_sig_bkg > cut) background += 1;
    if(i == ndata-1){
      sig_frac = (signal/float(ndata))*100;
      bkg_frac = (background/float(ndata))*100;
    }
    
    ///Getting Event Weights
    sig_event_weight = sig_neighbors/float(nsig);
    bkg_event_weight = bkg_neighbors/float(nbkg);
    event_weight     = sig_event_weight/(sig_event_weight + bkg_event_weight);
    
      
      ///Saving Current Results
      FME_out->Fill();
  }
  
  TFile *FME_Results = new TFile(Out_Name+"_Results.root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<"\n::::::::::::::::::::: Process Finished :::::::::::::::::::::"<<endl;
  cout<<"SigFrac = "<<sig_frac<<"%"<<endl;
  cout<<"BkgFrac = "<<bkg_frac<<"%"<<endl;

      
  seconds = difftime(stop,start);
  if(seconds < 60) elapsed_time = seconds;
  if(seconds >= 60){
    elapsed_time = seconds/60.;
    unity = "min.";
  }
  if(seconds >= 3600){
    elapsed_time = seconds/3600.;
    unity = "h";
  }
  cout<<" Time Consumed: "<<elapsed_time<<" "<<unity<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  
  return 0; //if process well finished
}




///:::::::::::::::::::::	FULL LEPTONIC + 2JETS	    :::::::::::::::::::::
int FME::FS_4l(TString Out_Name, TTree* Data_Tree, TTree* MC_Sig_Tree, TTree* MC_Bkg_Tree, TString Branch_Name){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
  
  cout<<"\n:::::::::::::::: Starting FastME Processing ::::::::::::::::"<<endl;
  cout<<"Model Chosen:  4l2j"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  Float_t Data[6][3][2], MC[6][3][2];
  Data_Tree->SetBranchAddress(Branch_Name,&Data);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Branch_Name,&MC);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Branch_Name,&MC);
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  double minDR_toSig, minDR_toBkg, prob_sig_bkg, sig_frac, event_weight;
  double bkg_frac, min_dr_sig, min_dr_bkg, sig_event_weight, bkg_event_weight;
  int signal = 0, background = 0, sig_neighbors, bkg_neighbors;
  TTree *FME_out = new TTree(Out_Name,"Fast Matrix Element Results");
  FME_out->SetDirectory(0);
  FME_out->Branch("minDR_toSig",&minDR_toSig);
  FME_out->Branch("minDR_toBkg",&minDR_toBkg);
  FME_out->Branch("SigFrac",&sig_frac);
  FME_out->Branch("BkgFrac",&bkg_frac);
  FME_out->Branch("P_SB",&prob_sig_bkg);
  FME_out->Branch("WSig_ToEvent",&sig_event_weight);
  FME_out->Branch("WBkg_ToEvent",&bkg_event_weight);
  FME_out->Branch("Event_Weight",&event_weight);

  
  for(int i=0; i<ndata; i++){
    if(i % (ndata/10) == 0) cout<<"Remaining Data: "<<ndata-i<<endl;
    Data_Tree->GetEntry(i);
    
    //Reseting Variables
    minDR_toSig   	= pedestal;
    minDR_toBkg   	= pedestal;
    prob_sig_bkg  	= pedestal;
    sig_frac      	= pedestal;
    bkg_frac	  	= pedestal;
    sig_event_weight 	= pedestal;
    bkg_event_weight 	= pedestal;
    min_dr_sig    	= 1.E15;
    min_dr_bkg    	= 1.E15;
    sig_neighbors 	= 0;
    bkg_neighbors 	= 0;
    
    for(int s=0; s<nsig; s++){
      MC_Sig_Tree->GetEntry(s);
      double dr_test = std_DR(6,Data,MC);
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
      if(dr_test < phs_radius) sig_neighbors += 1;
    }
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      double dr_test = std_DR(6,Data,MC);
      if(dr_test < min_dr_bkg) min_dr_bkg = dr_test;
      if(dr_test < phs_radius) bkg_neighbors += 1;
    }    
    
    ///Getting the discriminant value
    if(min_dr_sig == 0 && min_dr_bkg == 0){ 
      min_dr_sig = 1;
      min_dr_bkg = 1;
    }
      prob_sig_bkg = min_dr_sig/(min_dr_sig + min_dr_bkg);
      minDR_toSig = min_dr_sig;
      minDR_toBkg = min_dr_bkg;
      
    ///Checking amount of kind choices
    if(prob_sig_bkg < cut) signal += 1;
    if(prob_sig_bkg > cut) background += 1;
    if(i == ndata-1){
      sig_frac = (signal/float(ndata))*100;
      bkg_frac = (background/float(ndata))*100;
    }
    
    ///Getting Event Weights
    sig_event_weight = sig_neighbors/float(nsig);
    bkg_event_weight = bkg_neighbors/float(nbkg);
    event_weight     = sig_event_weight/(sig_event_weight + bkg_event_weight);
    
      
      ///Saving Current Results
      FME_out->Fill();
  }
  
  TFile *FME_Results = new TFile(Out_Name+"_Results.root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<"\n::::::::::::::::::::: Process Finished :::::::::::::::::::::"<<endl;
  cout<<"SigFrac = "<<sig_frac<<"%"<<endl;
  cout<<"BkgFrac = "<<bkg_frac<<"%"<<endl;

      
  seconds = difftime(stop,start);
  if(seconds < 60) elapsed_time = seconds;
  if(seconds >= 60){
    elapsed_time = seconds/60.;
    unity = "min.";
  }
  if(seconds >= 3600){
    elapsed_time = seconds/3600.;
    unity = "h";
  }
  cout<<" Time Consumed: "<<elapsed_time<<" "<<unity<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  
  return 0; //if process well finished
}




///:::::::::::::::::::::	SEMI-LEPTONIC CASE	   :::::::::::::::::::::
int FME::FS_lv2j(){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
    
  cout<<"\n:::::::::::::::: Starting FastME Processing ::::::::::::::::"<<endl;
  cout<<"Model Chosen:  lv2j"<<endl;
  cout<<"------------------------------------------------------------"<<endl;

  
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start=0;

    ///Electron
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    ///MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    ///Jet 1
    start = 2;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    ///Jet 2
    start = 3;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
  
  return 0; //if process well finished
}
///====================================================================================================