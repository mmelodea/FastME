//By Miquéias M. de Almeida 

#include <time.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>

#define l 6					//Number of particles (4l and 2jets)
#define c 4		  			//Number of variables (pT, eta, phi and E)
#define cut 0.5					//Threshold to separate events (ideal cut 0.5 - MC #Sig and #Bkg equals)

/****************************************************************************************************
* Calculate the distance between events using matrix with columms containing pT, eta, phi and lines *
* containing (l1-l4,j1,j2) the leptons indice. The leptons Data and MC are organized by closure     *
* constraining that 2 leptons came from the same vertex (constraining Ressonance)      		    *
****************************************************************************************************/
Float_t dr(Float_t Data[l][c], Float_t MC[l][c]){
  Float_t dpt2, deta2, dphi2, sum_dr1, sum_dr2;
  Float_t sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  do{
      sum_dr1 = 0;
      sum_dr2 = 0;
      
      dpt2    = pow(Data[start][0]-MC[start][0],2);
      deta2   = pow(Data[start][1]-MC[start][1],2);
      dphi2   = pow(Data[start][2]-MC[start][2],2);
      sum_dr1 = sqrt(dpt2 + deta2 + dphi2);
      
      dpt2    = pow(Data[start][0]-MC[start+1][0],2);
      deta2   = pow(Data[start][1]-MC[start+1][1],2);
      dphi2   = pow(Data[start][2]-MC[start+1][2],2);
      sum_dr2 = sqrt(dpt2 + deta2 + dphi2);
      
      if(sum_dr1 < sum_dr2){
	sum_dpt2  += pow(Data[start][0]-MC[start][0],2);
	sum_deta2 += pow(Data[start][1]-MC[start][1],2);
	sum_dphi2 += pow(Data[start][2]-MC[start][2],2);
	
	sum_dpt2  += pow(Data[start+1][0]-MC[start+1][0],2);
	sum_deta2 += pow(Data[start+1][1]-MC[start+1][1],2);
	sum_dphi2 += pow(Data[start+1][2]-MC[start+1][2],2);
      }
      
      if(sum_dr2 < sum_dr1){
	sum_dpt2  += pow(Data[start][0]-MC[start+1][0],2);
	sum_deta2 += pow(Data[start][1]-MC[start+1][1],2);
	sum_dphi2 += pow(Data[start][2]-MC[start+1][2],2);
	
	sum_dpt2  += pow(Data[start+1][0]-MC[start][0],2);
	sum_deta2 += pow(Data[start+1][1]-MC[start][1],2);
	sum_dphi2 += pow(Data[start+1][2]-MC[start][2],2);
      }

      start += 2; 
  }while(start < l);
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
//====================================================================================================


void FastME_Official(){
  
//---Timming process---
time_t date;
struct tm* time_info;
//---------------------
  
  int signal = 0, background = 0, inseparable = 0;
  Float_t data_stored[l][c], sig_stored[l][c], bkg_stored[l][c];
  Float_t dr_test, min_dr_sig, min_dr_bkg, prob_sig_bkg;
  
  TFile *data_file = TFile::Open("../Ntuples/Bkg.root");
  TTree *data = (TTree*)data_file->Get("VBF");
  data->SetBranchAddress("RECO_PARTICLE",&data_stored);
  Int_t ndata = (Int_t)data->GetEntries();
  
  TFile *sig_file = TFile::Open("../Ntuples/Sig2.root");
  TTree *sig = (TTree*)sig_file->Get("VBF");
  sig->SetBranchAddress("RECO_PARTICLE",&sig_stored);
  Int_t nsig = (Int_t)sig->GetEntries();
  
  TFile *bkg_file = TFile::Open("../Ntuples/Bkg_33e4.root");
  TTree *bkg = (TTree*)bkg_file->Get("VBF");
  bkg->SetBranchAddress("RECO_PARTICLE",&bkg_stored);
  Int_t nbkg = (Int_t)bkg->GetEntries();
  
  TTree *discriminator = new TTree("discriminator","Results from Discriminator");
  discriminator->SetDirectory(0);
  discriminator->Branch("prob_sig_bkg",&prob_sig_bkg,"prob_sig_bkg/F");  
  
  TH1D *dr_sig = new TH1D("dr_sig","",100,0,500);
  TH1D *dr_bkg = new TH1D("dr_bkg","",100,0,500);
  TH2D *prob_sig_bkg_Vl = new TH2D("prob_sig_bkg_Vl","",125,0,250,125,0,250);
  TH2D *prob_sig_bkg_Ev = new TH2D("prob_sig_bkg_Ev","",125,0,250,125,0,250);
  TH1D *sig_frac = new TH1D("sig_frac","",10000,0,101);
  TH1D *bkg_frac = new TH1D("bkg_frac","",10000,0,101);
  TH1D *insep_frac = new TH1D("insep_frac","",10000,0,101);

  //Controls to number of events
  ndata = 5.E3;
  nsig  = 5.E3;
  nbkg  = 5.E3;
      
  //Print on PC screen information about finish of processes
  time(&date);
  time_info = localtime(&date);
  printf("%s\n","_________________________________________________________");
  printf("%s\n","                Starting FastME Analysis                 ");
  printf("\t\t%s",asctime(time_info));
  printf("%s\n","Remind you, different number of MC samples reduce the purity!");
  printf("%s\n","_________________________________________________________");
  printf("%s\t %i\n","nData =",ndata);
  printf("%s\t %i\n","nSig  =",nsig);
  printf("%s\t %i\n","nBkg  =",nbkg);
  printf("%s\n","---------------------------------------------------------");

  for(Int_t i=0; i<ndata; i++){
    if(i % (ndata/10) == 0) cout<<"Remaining Data: "<<ndata-i<<endl;
    data->GetEntry(i);
    
    //Searching the minimum distance between Signal and Data event
    min_dr_sig = 1.E15;
    for(Int_t j=0; j<nsig; j++){
      sig->GetEntry(j);
      dr_test = dr(data_stored,sig_stored);
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
    }
      
    //Searching the minimum distance between Background and Data event
    min_dr_bkg = 1.E15;
    for(Int_t j=0; j<nbkg; j++){
      bkg->GetEntry(j);
      dr_test = dr(data_stored,bkg_stored);
      if(dr_test < min_dr_bkg) min_dr_bkg = dr_test;
    }
    
    /*********************************************************
    * Checking the probabilty to be Signal or Background     *
    * if prob_sig_bkg = 0.5 event can't be classificated     *
    * if prob_sig_bkg < 0.5 event is signal		     *
    * if prob_sig_bkg > 0.5 event is background		     *
    * ****************************************************** *
    * The inseparability occurs when min_dr_sig = min_dr_bkg *
    * and in this case pSB should be 0.5, otherwise if event *
    * is signal, min_dr_bkg is large than min_dr_sig and pSB *
    * can't be large then 0.5, and the same way to bkg event * 
    *********************************************************/
    if(min_dr_sig == 0 && min_dr_sig == 0){ 
      min_dr_sig = 1;
      min_dr_bkg = 1;
    }
    prob_sig_bkg = min_dr_sig/(min_dr_sig + min_dr_bkg); 
         
    //Checking amount of correct/wrong choices
    if(prob_sig_bkg < cut) signal += 1;
    if(prob_sig_bkg > cut) background += 1;
    if(prob_sig_bkg == cut) inseparable += 1;

    //Saving results
    dr_sig->Fill(min_dr_sig);
    dr_bkg->Fill(min_dr_bkg);
    prob_sig_bkg_Vl->Fill(min_dr_sig,min_dr_bkg,prob_sig_bkg);
    prob_sig_bkg_Ev->Fill(min_dr_sig,min_dr_bkg,1);
    discriminator->Fill();
  }
  
  //Saving the amount of event types indentificated
  float fsig = (signal*100)/float(ndata);
  float fbkg = (background*100)/float(ndata);
  float fins = (inseparable*100)/float(ndata);
  sig_frac->Fill(fsig);
  bkg_frac->Fill(fbkg);
  insep_frac->Fill(fins);
  
  TFile *file = TFile::Open("Teste_Bkg.root","recreate");
  discriminator->Write();
  dr_sig->Write();
  dr_bkg->Write();
  prob_sig_bkg_Vl->Write();
  prob_sig_bkg_Ev->Write();
  sig_frac->Write();
  bkg_frac->Write();
  insep_frac->Write();
  file->Close();
  
  //Print on PC screen information about finish of processes
  time(&date);
  time_info = localtime(&date);
  printf("%s\n","_________________________________________________________");
  printf("%s\n","                 FastME Analysis Finished                ");
  printf("%s\n","---------------------------------------------------------");
  printf("%s\t %.1f%s\n","Signal      =",fsig,"%");
  printf("%s\t %.1f%s\n","Background  =",fbkg,"%");
  printf("%s\t %.1f%s\n","Inseparable =",fins,"%");
  printf("%s\n","---------------------------------------------------------");
  printf("\t\t%s",asctime(time_info));
  printf("%s\n","_________________________________________________________");
  
}//End main program
