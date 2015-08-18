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
  Int_t ndata = (Int_t)data->GetEntries();
  
  TFile *sig_file = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_sig_formated_full.root");
  TTree *sig = (TTree*)sig_file->Get("VBF");
  sig->SetBranchAddress("RECO_PARTICLE",&sig_stored);
  Int_t nsig = (Int_t)sig->GetEntries();

  TFile *bkg1_file = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_bkg_qqzz_part1.root");
  TTree *bkg1 = (TTree*)bkg1_file->Get("VBF");
  bkg1->SetBranchAddress("RECO_PARTICLE",&bkg1_stored);
  Int_t nbkg1 = (Int_t)bkg1->GetEntries();  
    
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
  //printf("%s\t %i\n","nBkg2 =",nbkg2);
  //printf("%s\t %.3f\n","Cut   =",cut);
  printf("%s\n","---------------------------------------------------------");

  FME *teste = new FME();
  teste->launch_FME("4l","FastME",data,sig,bkg1,"RECO_PARTICLE");
  
  //Print on PC screen information about finish of processes
  time(&date);
  time_info = localtime(&date);
  printf("%s\n","_________________________________________________________");
  printf("%s\n","                 FastME Analysis Finished                ");
  printf("%s\n","---------------------------------------------------------");
  printf("\t\t%s",asctime(time_info));
  printf("%s\n","_________________________________________________________");
  
}//End main program