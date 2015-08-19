///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	FAST MATRIX ELEMENT CLASS MEMBERS DEFINITIONS	::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::          Author: Miquéias M. de Almeida 		::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#include "../PhsDrComputers/FS4l_DrComputers.cxx"
#include "../PhsDrComputers/FS4l2j_DrComputers.cxx"
#include "../PhsDrComputers/FSlv2j_DrComputers.cxx"
#include "FastME.h"
#include <iostream>
#include <string>
#include <ctime>
#include <exception>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#define pedestal 	-99					///Reset Value to Variables
#define cut 	 	0.5					///Threshold to Separate Events (Ideal Cut 0.5 - MC #Sig and #Bkg Equals)

using namespace std;


int FME::launchFME(TString Final_State, TString Model, TString Out_Name, TString Data_Path, TString MC_Sig_Path, 
		   TString MC_Bkg_Path, TString Tree_Name, TString Branch_Name, TString Resolution, TString PHS_Radius){
  const int nFS = 3;
  int FS = -1;
  TString FState[nFS] = {"4l","4l2j","lv2j"};
       if(Final_State == FState[0]) FS = 0;
  else if(Final_State == FState[1]) FS = 1;
  else if(Final_State == FState[2]) FS = 2;
  else{
    cout<<"[Error] Final State '"<<Final_State<<"' not defined!"<<endl;
    cout<<"Possible candidates are: ";
    for(int nM=0; nM<nFS; nM++) cout<<FState[nM]<<", ";
    cout<<endl;
    throw exception();
  }
  
  float PHS_radius = PHS_Radius.Atof();
  if(PHS_radius < 0){
    cout<<"[Error] Phase Space Radius threshold is negative!"<<endl;
    throw exception();
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
      return FS_4l(Model, Out_Name, Data_Tree, MC_Sig_Tree, MC_Bkg_Tree, Branch_Name, Resolution, PHS_radius);
      break;
      
    case 1:
      return FS_4l2j(Model, Out_Name, Data_Tree, MC_Sig_Tree, MC_Bkg_Tree, Branch_Name, Resolution, PHS_radius);
      break;
      
    case 2:
      return FS_lv2j(Model, Out_Name, Data_Tree, MC_Sig_Tree, MC_Bkg_Tree, Branch_Name, Resolution, PHS_radius);
      break;
    
    default:
      return FS_4l(Model, Out_Name, Data_Tree, MC_Sig_Tree, MC_Bkg_Tree, Branch_Name, Resolution, PHS_radius);
      break;
  }
  return 0;
}


///:::::::::::::::::::::	FULL LEPTONIC CASE	   :::::::::::::::::::::
int FME::FS_4l(TString Model, TString Out_Name, TTree* Data_Tree, TTree* MC_Sig_Tree, TTree* MC_Bkg_Tree,
	       TString Branch_Name, TString Resolution, Double_t PHS_radius){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
  
  cout<<"\n:::::: Starting FastME Processing ::::::"<<endl;
  cout<<":: Final State:    "<<"4l"<<endl;
  cout<<":: Model Chosen:   "<<Model<<endl;
  cout<<":: PHS Radius:     "<<PHS_radius<<endl;
  cout<<":: Use Resolution: "<<Resolution<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  Float_t Data[4][3][2], MC[4][3][2];
  Data_Tree->SetBranchAddress(Branch_Name,&Data);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Branch_Name,&MC);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Branch_Name,&MC);
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  double dr_test, minDR_toSig, minDR_toBkg, prob_sig_bkg, sig_frac, event_weight;
  double bkg_frac, min_dr_sig, min_dr_bkg, sig_event_weight, bkg_event_weight;
  int signal = 0, background = 0, sig_neighbors, bkg_neighbors;
  TTree *FME_out = new TTree("FastME_Results","Fast Matrix Element Results");
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
    if(i % (ndata/10) == 0) cout<<":: Remaining Data: "<<ndata-i<<endl;
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
      
            if(Model == "DR_Order_Res")   dr_test = FS4l_Res_DrOrder(Data,MC,Resolution);
       else if(Model == "DR_Order_noRes") dr_test = FS4l_noRes_DrOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_Res")   dr_test = FS4l_Res_PtOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_noRes") dr_test = FS4l_noRes_PtOrder(Data,MC,Resolution);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
      if(dr_test < PHS_radius) sig_neighbors += 1;
    }
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      
            if(Model == "DR_Order_Res")   dr_test = FS4l_Res_DrOrder(Data,MC,Resolution);
       else if(Model == "DR_Order_noRes") dr_test = FS4l_noRes_DrOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_Res")   dr_test = FS4l_Res_PtOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_noRes") dr_test = FS4l_noRes_PtOrder(Data,MC,Resolution);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      
      if(dr_test < min_dr_bkg) min_dr_bkg = dr_test;
      if(dr_test < PHS_radius) bkg_neighbors += 1;
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
  
  TString path = "../interface/";
  TFile *FME_Results = new TFile(path+Out_Name+".root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<"::::::::::: Process Finished :::::::::::"<<endl;
  cout<<":: SigFrac =      "<<sig_frac<<"%"<<endl;
  cout<<":: BkgFrac =      "<<bkg_frac<<"%"<<endl;

      
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
  cout<<":: Time Consumed: "<<elapsed_time<<unity<<endl;
  cout<<"::--------------------------------------"<<endl;
  
  
  return 0; //if process well finished
}




///:::::::::::::::::::::	FULL LEPTONIC + 2JETS	    :::::::::::::::::::::
int FME::FS_4l2j(TString Model, TString Out_Name, TTree* Data_Tree, TTree* MC_Sig_Tree, TTree* MC_Bkg_Tree,
		 TString Branch_Name, TString Resolution, Double_t PHS_radius){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
  
  cout<<"\n:::::::::::::::: Starting FastME Processing ::::::::::::::::"<<endl;
  cout<<":: Final State:    "<<"4l2j"<<endl;
  cout<<":: Model Chosen:   "<<Model<<endl;
  cout<<":: PHS Radius:     "<<PHS_radius<<endl;
  cout<<":: Use Resolution: "<<Resolution<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  Float_t Data[6][3][2], MC[6][3][2];
  Data_Tree->SetBranchAddress(Branch_Name,&Data);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Branch_Name,&MC);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Branch_Name,&MC);
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  double dr_test, minDR_toSig, minDR_toBkg, prob_sig_bkg, sig_frac, event_weight;
  double bkg_frac, min_dr_sig, min_dr_bkg, sig_event_weight, bkg_event_weight;
  int signal = 0, background = 0, sig_neighbors, bkg_neighbors;
  TTree *FME_out = new TTree("FastME_Results","Fast Matrix Element Results");
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
      
            if(Model == "DR_Order_Res")   dr_test = FS4l2j_Res_DrOrder(Data,MC,Resolution);
       else if(Model == "DR_Order_noRes") dr_test = FS4l2j_noRes_DrOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_Res")   dr_test = FS4l2j_Res_PtOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_noRes") dr_test = FS4l2j_noRes_PtOrder(Data,MC,Resolution);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
      if(dr_test < PHS_radius) sig_neighbors += 1;
    }
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      
            if(Model == "DR_Order_Res")   dr_test = FS4l2j_Res_DrOrder(Data,MC,Resolution);
       else if(Model == "DR_Order_noRes") dr_test = FS4l2j_noRes_DrOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_Res")   dr_test = FS4l2j_Res_PtOrder(Data,MC,Resolution);
       else if(Model == "PT_Order_noRes") dr_test = FS4l2j_noRes_PtOrder(Data,MC,Resolution);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      
      if(dr_test < min_dr_bkg) min_dr_bkg = dr_test;
      if(dr_test < PHS_radius) bkg_neighbors += 1;
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
  
  TString path = "../interface/";
  TFile *FME_Results = new TFile(path+Out_Name+".root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<"::::::::::::::::::::: Process Finished :::::::::::::::::::::"<<endl;
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
int FME::FS_lv2j(TString Model, TString Out_Name, TTree* Data_Tree, TTree* MC_Sig_Tree, TTree* MC_Bkg_Tree,
		 TString Branch_Name, TString Resolution, Double_t PHS_radius){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
    
  cout<<"\n:::::::::::::::: Starting FastME Processing ::::::::::::::::"<<endl;
  cout<<":: Final State:    "<<"lv2j"<<endl;
  cout<<":: Model Chosen:   "<<Model<<endl;
  cout<<":: PHS Radius:     "<<PHS_radius<<endl;
  cout<<":: Use Resolution: "<<Resolution<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  Float_t Data[4][3][2], MC[4][3][2];
  Data_Tree->SetBranchAddress(Branch_Name,&Data);
  int ndata = Data_Tree->GetEntries();
  MC_Sig_Tree->SetBranchAddress(Branch_Name,&MC);
  int nsig = MC_Sig_Tree->GetEntries();
  MC_Bkg_Tree->SetBranchAddress(Branch_Name,&MC);
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  double dr_test, minDR_toSig, minDR_toBkg, prob_sig_bkg, sig_frac, event_weight;
  double bkg_frac, min_dr_sig, min_dr_bkg, sig_event_weight, bkg_event_weight;
  int signal = 0, background = 0, sig_neighbors, bkg_neighbors;
  TTree *FME_out = new TTree("FastME_Results","Fast Matrix Element Results");
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
      
            if(Model == "Jet_DR_Order") dr_test = FSlv2j_Jet_DrOrder(Data,MC,Resolution);
       else if(Model == "Jet_PT_Order") dr_test = FSlv2j_Jet_PtOrder(Data,MC,Resolution);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
      if(dr_test < PHS_radius) sig_neighbors += 1;
    }
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      
            if(Model == "Jet_DR_Order") dr_test = FSlv2j_Jet_DrOrder(Data,MC,Resolution);
       else if(Model == "Jet_PT_Order") dr_test = FSlv2j_Jet_PtOrder(Data,MC,Resolution);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      
      if(dr_test < min_dr_bkg) min_dr_bkg = dr_test;
      if(dr_test < PHS_radius) bkg_neighbors += 1;
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
  
  TString path = "../interface/";
  TFile *FME_Results = new TFile(path+Out_Name+".root","recreate");
  FME_out->Write();
  FME_Results->Close();
  
  time(&stop);
  cout<<"::::::::::::::::::::: Process Finished :::::::::::::::::::::"<<endl;
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
///====================================================================================================