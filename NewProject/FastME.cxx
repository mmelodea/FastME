///FastME Constructor
///Author: Miquéias M. de Almeida 

#include "FastME.h"
#include <iostream>
#include <ctime>
#include <exception>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#define pedestal -99					///Reset Value to Variables
#define cut 	 0.5					///Threshold to Separate Events (Ideal Cut 0.5 - MC #Sig and #Bkg Equals)
#define phs_radius 300

using namespace std;


///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::: PHASE ESPACE DATA-MC DISTANCE COMPUTER FOR EQUAL RESSONANCE COMPONENTS :::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

///DISTANCE ORDERING + RESSONANCE CRITERY
double std_DR(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  do{
      sum_dr1 = 0;
      sum_dr2 = 0;
      
      fdpt2    = pow( (Data[start][0][0]-MC[start][0][0])/Data[start][0][1] ,2 );
      fdeta2   = pow( (Data[start][1][0]-MC[start][1][0])/Data[start][1][1] ,2 );
      fdphi2   = pow( (Data[start][2][0]-MC[start][2][0])/Data[start][2][1] ,2 );
      sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
      sdpt2    = pow( (Data[start][0][0]-MC[start+1][0][0])/Data[start][0][1] ,2 );
      sdeta2   = pow( (Data[start][1][0]-MC[start+1][1][0])/Data[start][1][1] ,2 );
      sdphi2   = pow( (Data[start][2][0]-MC[start+1][2][0])/Data[start][2][1] ,2 );
      sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
    if(sum_dr1 < sum_dr2){
      sum_dpt2  += fdpt2;
      sum_deta2 += fdeta2;
      sum_dphi2 += fdphi2;
  
      sum_dpt2  += pow( (Data[start+1][0][0]-MC[start+1][0][0])/Data[start+1][0][1] ,2 );
      sum_deta2 += pow( (Data[start+1][1][0]-MC[start+1][1][0])/Data[start+1][1][1] ,2 );
      sum_dphi2 += pow( (Data[start+1][2][0]-MC[start+1][2][0])/Data[start+1][2][1] ,2 );
    }
      
    if(sum_dr2 < sum_dr1){
      sum_dpt2  += sdpt2;
      sum_deta2 += sdeta2;
      sum_dphi2 += sdphi2;

      sum_dpt2  += pow( (Data[start+1][0][0]-MC[start][0][0])/Data[start+1][0][1] ,2 );
      sum_deta2 += pow( (Data[start+1][1][0]-MC[start][1][0])/Data[start+1][1][1] ,2 );
      sum_dphi2 += pow( (Data[start+1][2][0]-MC[start][2][0])/Data[start+1][2][1] ,2 );
    }

    start += 2; 
  }while(start < 4);
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  
  return event_distance;
}


///DISTANCE ORDERING + NO RESSONANCE CRITERY
double noRes_DR(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double dpt2, deta2, dphi2, sum_dr = 0, limiar;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int chosen = -1, p[] = {-1,-1,-1,-1,-1,-1};
  double vchosen[] = {-1,-1,-1};
  
  //Choose the l/j1-Data most close to l/j1-MC
  for(int i=0; i<4; i++){
    limiar = 1.E15;
    for(int j=0; j<4; j++){
      //if((j > 3 && i < 4) || (j < 4 && i > 3)) continue;       		                ///Avoid leptons/jets mixture
      if(j == p[0] || j == p[1] || j == p[2]) continue;	///Avoid recounting
      //if(j == p[0] || j == p[1] || j == p[2] || j == p[3] || j == p[4]) continue;	///Avoid recounting
      
      dpt2  = pow( (Data[i][0][0]-MC[j][0][0])/Data[i][0][1] ,2 );
      deta2 = pow( (Data[i][1][0]-MC[j][1][0])/Data[i][1][1] ,2 );
      dphi2 = pow( (Data[i][2][0]-MC[j][2][0])/Data[i][2][1] ,2 );
      sum_dr = sqrt(dpt2 + deta2 + dphi2);
      
      if(sum_dr < limiar){
	limiar     = sum_dr;
	chosen     = j;
	vchosen[0] = dpt2;
	vchosen[1] = deta2;
	vchosen[2] = dphi2;
      }
    }
    p[i] = chosen;
    
    //Makes the sum of deltapT2, deta2 and dphi2 to get the events distance
    sum_dpt2  += vchosen[0];
    sum_deta2 += vchosen[1];
    sum_dphi2 += vchosen[2];
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}


///PT ORDERING + RESSONANCE CRITERY
double pt_Res_DR(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double data_limiar, mc_limiar;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int data_ch = -1, mc_ch = -1, data_p[] = {-1,-1,-1,-1,-1,-1}, mc_p[] = {-1,-1,-1,-1,-1,-1};
  
  //Organizing leptons by pT
  for(int i=0; i<4; i++){
    data_limiar = 0.;
    mc_limiar   = 0.;
        
    for(int j=0; j<4; j++){
      //if((j > 3 && i < 4) || (j < 4 && i > 3)) continue;       		        ///Avoid leptons/jets mixture
      
      if(j != data_p[0] && j != data_p[1] && j != data_p[2] && j != data_p[4]){         ///MC avoid recounting
	if(Data[j][0][0] > data_limiar){
	  data_limiar = Data[j][0][0];
	  data_ch = j;
	}
      }
      if(j != mc_p[0] && j != mc_p[1] && j != mc_p[2] && j != mc_p[4]){			///Data avoid recounting
	if(MC[j][0][0] > mc_limiar){
	  mc_limiar = MC[j][0][0];
	  mc_ch = j;
	}
      }
    }
    data_p[i] = data_ch;
    mc_p[i] = mc_ch;

    //Makes the sum over all final state particles
    sum_dpt2  += pow( (Data[data_ch][0][0]-MC[mc_ch][0][0])/Data[data_ch][0][1] ,2 );
    sum_deta2 += pow( (Data[data_ch][1][0]-MC[mc_ch][1][0])/Data[data_ch][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_ch][2][0]-MC[mc_ch][2][0])/Data[data_ch][2][1] ,2 );
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}


///PT ORDERING + NO RESSONANCE CRITERY
double pt_noRes_DR(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  Float_t sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int data_st = -1, data_nd = -1, mc_st = -1, mc_nd = -1, start = 0;
 
  //Organizing leptons by pT
  for(int i=0; i<2; i++){
    
    //Organizing data
    if(Data[start][0][0] > Data[start+1][0][0]){ data_st = start; data_nd = start+1; }
    if(Data[start][0][0] < Data[start+1][0][0]){ data_st = start+1; data_nd = start; }
      
    //Organizing MC
    if(MC[start][0][0] > MC[start+1][0][0]){ mc_st = start; mc_nd = start+1; }
    if(MC[start][0][0] < MC[start+1][0][0]){ mc_st = start+1; mc_nd = start; }

    //Makes the sum over all final state particles
    sum_dpt2  += pow( (Data[data_st][0][0]-MC[mc_st][0][0])/Data[data_st][0][1] ,2 );
    sum_deta2 += pow( (Data[data_st][1][0]-MC[mc_st][1][0])/Data[data_st][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_st][2][0]-MC[mc_st][2][0])/Data[data_st][2][1] ,2 );
    
    sum_dpt2  += pow( (Data[data_nd][0][0]-MC[mc_nd][0][0])/Data[data_nd][0][1] ,2 );
    sum_deta2 += pow( (Data[data_nd][1][0]-MC[mc_nd][1][0])/Data[data_nd][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_nd][2][0]-MC[mc_nd][2][0])/Data[data_nd][2][1] ,2 );
    
    start += 2;
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///==========================================================================================




///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	FAST MATRIX ELEMENT CLASS MEMBERS DEFINITIONS	::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int FME::launchFME(TString Final_State, TString Model, TString Out_Name, TString Data_Path, TString MC_Sig_Path, TString MC_Bkg_Path, TString Tree_Name, TString Branch_Name){
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
      return FS_4l(Model,Out_Name,Data_Tree,MC_Sig_Tree,MC_Bkg_Tree,Branch_Name);
      break;
/*    
    case 1:
      return FS_4l2j(Model,Out_Name,Data_Tree,MC_Sig_Tree,MC_Bkg_Tree,Branch_Name);
      break;
      
    case 2:
      return FS_lv2j_DR(Model,TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
      break;
    
    default:
      return FS_4l(Model,TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
      break;
      */
  }
  return 0;
}


///:::::::::::::::::::::	FULL LEPTONIC CASE	   :::::::::::::::::::::
int FME::FS_4l(TString Model, TString Out_Name, TTree* Data_Tree, TTree* MC_Sig_Tree, TTree* MC_Bkg_Tree, TString Branch_Name){
  ///For timming the process
  time_t start, stop;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///-----------------------
  
  cout<<"\n:::::: Starting FastME Processing ::::::"<<endl;
  cout<<":: Final State:   "<<"4l"<<endl;
  cout<<":: Model Chosen:  "<<Model<<endl;
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
      
            if(Model == "DR_Order_Res")   dr_test = std_DR(Data,MC);
       else if(Model == "DR_Order_noRes") dr_test = noRes_DR(Data,MC);
       else if(Model == "PT_Order_Res")   dr_test = pt_Res_DR(Data,MC);
       else if(Model == "PT_Order_noRes") dr_test = pt_noRes_DR(Data,MC);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
      if(dr_test < phs_radius) sig_neighbors += 1;
    }
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      
            if(Model == "DR_Order_Res")   dr_test = std_DR(Data,MC);
       else if(Model == "DR_Order_noRes") dr_test = noRes_DR(Data,MC);
       else if(Model == "PT_Order_Res")   dr_test = pt_Res_DR(Data,MC);
       else if(Model == "PT_Order_noRes") dr_test = pt_noRes_DR(Data,MC);
       else{
	 cout<<"Model '"<<Model<<"' is not defined!"<<endl;
	 throw exception();
      }
      
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
  cout<<"\n::::::::::: Process Finished :::::::::::"<<endl;
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


/*

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
*/
///====================================================================================================