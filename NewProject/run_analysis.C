#include "FastME.cxx"
#include "FastME.h"

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

#define cut 0.5					//Controls the Sig/Bkg discrimination (0.5: ideal case)

void run_analysis(){
  
//--Timming processes--
time_t date;
struct tm* time_info;
//---------------------
  
  int signal = 0, background = 0, inseparable = 0;
  Float_t *data_stored, *sig_stored, *bkg1_stored, *bkg2_stored, reco, weight, mela, kd;
  Float_t dr_test, min_dr_sig, min_dr_bkg1, min_dr_bkg2, min_dr_bkg, prob_sig_bkg, dr_sig, dr_bkg;
  
  TFile *data_file = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_data_formated.root");
  TTree *data = (TTree*)data_file->Get("VBF");
  data->SetBranchAddress("RECO_PARTICLE",&data_stored);
  //data->SetBranchAddress("reco",&reco);
  //data->SetBranchAddress("weight",&weight);
  //data->SetBranchAddress("f_KD",&mela);
  Int_t ndata = (Int_t)data->GetEntries();
  
  TFile *sig_file = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_sig_formated_full.root");
  TTree *sig = (TTree*)sig_file->Get("VBF");
  sig->SetBranchAddress("RECO_PARTICLE",&sig_stored);
  Int_t nsig = (Int_t)sig->GetEntries();

  TFile *bkg1_file = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_bkg_qqzz_part1.root");
  TTree *bkg1 = (TTree*)bkg1_file->Get("VBF");
  bkg1->SetBranchAddress("RECO_PARTICLE",&bkg1_stored);
  Int_t nbkg1 = (Int_t)bkg1->GetEntries();
  
  TFile *bkg2_file = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_bkg_ggzz_part1.root");
  TTree *bkg2 = (TTree*)bkg2_file->Get("VBF");
  bkg2->SetBranchAddress("RECO_PARTICLE",&bkg2_stored);
  Int_t nbkg2 = (Int_t)bkg2->GetEntries();
  
  TH1D *frac_sig = new TH1D("frac_sig","",12000,-1000,11000);
  TH1D *frac_bkg = new TH1D("frac_bkg","",12000,-1000,11000);

  TTree *discriminator = new TTree("discriminator","Results from Discriminator");
  discriminator->SetDirectory(0);
  discriminator->Branch("prob_sig_bkg",&prob_sig_bkg,"prob_sig_bkg/F");
  discriminator->Branch("kd",&kd,"kd/F");
  discriminator->Branch("reco",&reco,"reco/F");
  discriminator->Branch("dr_sig",&dr_sig,"dr_sig/F");
  discriminator->Branch("dr_bkg",&dr_bkg,"dr_bkg/F");
  
  
  TH2D *prob_sig_bkg_Vl = new TH2D("prob_sig_bkg_Vl","",125,0,250,125,0,250);
  TH2D *prob_sig_bkg_Ev = new TH2D("prob_sig_bkg_Ev","",125,0,250,125,0,250);
  TH1D *sig_frac = new TH1D("sig_frac","",10000,0,101);
  TH1D *bkg_frac = new TH1D("bkg_frac","",10000,0,101);
  TH1D *insep_frac = new TH1D("insep_frac","",10000,0,101);
  TH1D *reco_sig = new TH1D("reco_sig","",70,100,800);
  TH1D *reco_bkg = new TH1D("reco_bkg","",70,100,800);
  TH1D *orig_sig = new TH1D("orig_sig","",70,100,800);
  TH1D *pSB = new TH1D("pSB","",50,0,1);
  TH1D *KD = new TH1D("KD","",50,0,1);
  TH2D *m4l_pSB = new TH2D("m4l_pSB","",140,100,800,50,0,1);
  TH2D *m4l_KD = new TH2D("m4l_KD","",140,100,800,50,0,1);

  //Controls the number of events for test
  //ndata = 1.E4;
  nsig  = nbkg1+nbkg2;
  //nbkg1 = nsig/2;//int((nbkg1*nsig)/float(nbkg1+nbkg2));
  //nbkg2 = nsig/2;//int((nbkg2/float(nbkg1+nbkg2))*nsig);
  //Float_t cut = cut_cor(nsig,nbkg);
    
  //Print on PC screen information about finish of processes
  time(&date);
  time_info = localtime(&date);
  printf("%s\n","_________________________________________________________");
  printf("%s\n","                 Starting FastME Analysis                ");
  printf("\t\t%s",asctime(time_info));
  printf("%s\n","_________________________________________________________");
  printf("%s\t %i\n","nData =",ndata);
  printf("%s\t %i\n","nSig  =",nsig);
  printf("%s\t %i\n","nBkg1 =",nbkg1);
  printf("%s\t %i\n","nBkg2 =",nbkg2);
  //printf("%s\t %.3f\n","Cut   =",cut);
  printf("%s\n","---------------------------------------------------------");

  FME *teste = new FME();
  for(Int_t i=0; i<ndata; i++){
    if(i % (ndata/10) == 0) cout<<"Remaining Data: "<<ndata-i<<endl;
    data->GetEntry(i);
    
    //Checking medium distance between Signal and Data events
    min_dr_sig = 1.E15;
    for(Int_t j=0; j<nsig; j++){
      sig->GetEntry(j);
      dr_test = teste->compute_DR("4l",&data_stored,&sig_stored);
      if(dr_test < min_dr_sig) min_dr_sig = dr_test;
    }
      
    //Checking medium distance between Background and Data events
    min_dr_bkg1 = 1.E15;
    for(Int_t j=0; j<nbkg1; j++){
      bkg1->GetEntry(j);
      dr_test = teste->compute_DR("4l",&data_stored,&bkg1_stored);
      if(dr_test < min_dr_bkg1) min_dr_bkg1 = dr_test;
    }
    
    //Checking medium distance between Background and Data events
    min_dr_bkg2 = 1.E15;
    for(Int_t j=0; j<nbkg2; j++){
      bkg2->GetEntry(j);
      dr_test = teste->compute_DR("4l",&data_stored,&bkg2_stored);
      if(dr_test < min_dr_bkg2) min_dr_bkg2 = dr_test;
    }
    min_dr_bkg = (min_dr_bkg1 < min_dr_bkg2)? min_dr_bkg1:min_dr_bkg2;
    /*********************************************************
    *   Checking the probabilty to be Signal or Background   *
    *   if prob_sig_bkg = 0.5 event can't be classificated   *
    *   if prob_sig_bkg < 0.5 event is signal		     *
    *   if prob_sig_bkg > 0.5 event is background	     *
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
    prob_sig_bkg = min_dr_sig/(min_dr_sig+min_dr_bkg);
    pSB->Fill(prob_sig_bkg);
         
    //Checking amount of correct/wrong choices
    if(prob_sig_bkg < cut) signal += 1;
    if(prob_sig_bkg > cut) background += 1;
    if(prob_sig_bkg == cut) inseparable += 1;
    if(prob_sig_bkg < cut) reco_sig->Fill(reco);
    if(prob_sig_bkg > cut) reco_bkg->Fill(reco);

    //Saving results
    dr_sig = min_dr_sig;
    dr_bkg = min_dr_bkg;
    prob_sig_bkg_Vl->Fill(min_dr_sig,min_dr_bkg,prob_sig_bkg);
    prob_sig_bkg_Ev->Fill(min_dr_sig,min_dr_bkg,1);
    m4l_pSB->Fill(reco,prob_sig_bkg);
    m4l_KD->Fill(reco,mela);
    KD->Fill(mela);
    kd = mela;
    discriminator->Fill();
    orig_sig->Fill(reco);
  }
  
  //Getting the amount correct/wrong choices
  float fsig = (signal*100)/float(ndata);
  float fbkg = (background*100)/float(ndata);
  float fins = (inseparable*100)/float(ndata);
  //frac_sig->Fill(h,fsig);
  //frac_bkg->Fill(h,fbkg);
  
  //cout<<"fsig: "<<fsig<<"   fbkg: "<<fbkg<<"   Remaining: "<<10001-h<<endl;
  //if(h<=300) h += 20;
  //if(h>300) h += 1000;
//}while(h<10001);
  sig_frac->Fill(fsig);
  bkg_frac->Fill(fbkg);
  insep_frac->Fill(fins);
  
  TFile *file = TFile::Open("first_teste.root","recreate");
  //frac_sig->Write();
  //frac_bkg->Write();
  discriminator->Write();
  prob_sig_bkg_Vl->Write();
  prob_sig_bkg_Ev->Write();
  sig_frac->Write();
  bkg_frac->Write();
  insep_frac->Write();
  reco_sig->Write();
  reco_bkg->Write();
  orig_sig->Write();
  pSB->Write();
  KD->Write();
  m4l_pSB->Write();
  m4l_KD->Write();
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