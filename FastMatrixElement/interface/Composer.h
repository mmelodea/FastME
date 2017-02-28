///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::[ Composer -  Event generator based on min DR ]::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Composer_h
#define Composer_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <cmath>
#include <iomanip>

#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TRandom3.h>

///Headers to TProcPool
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>



Int_t get_number_of_entries(TString file_name, TString tree_name){
  TFile *file = TFile::Open(file_name);
  TTree *tree = (TTree*)file->Get(tree_name);
  Int_t nentries = tree->GetEntries();
  tree->Delete();
  file->Close();
  
  return nentries;
}


double get_event_distance(int nMcParticles,
			  std::vector<Double_t>& Pt_1, std::vector<Double_t>& Eta_1, std::vector<Double_t>& Phi_1,
			  std::vector<Double_t>* Pt_2, std::vector<Double_t>* Eta_2, std::vector<Double_t>* Phi_2){
  
  //Loop over the particles in the data event
  std::vector<int> McFlag(nMcParticles,0);
  Double_t Min_dPt2_dEta2_dPhi2 = 0, Sum_dPt2_dEta2_dPhi2 = 0;
  for(int imc1 = 0; imc1 < nMcParticles; ++imc1){
    Double_t min_particles_distance = 1.e15;
    Int_t sel_mc_part = -1;
    
    //Loop over the particles in the MC event
    for(Int_t imc2 = 0; imc2 < nMcParticles; ++imc2){
      ///Avoid MC flaged as "used"
      if(McFlag[imc2] == 1) continue;

      ///Compute preliminary particles distance
      Double_t dPt2  = pow( (Pt_1[imc1]-(*Pt_2)[imc2]), 2 );	      
      Double_t dEta2 = pow( (Eta_1[imc1]-(*Eta_2)[imc2]), 2 );
      Double_t dPhi2 = pow( (Phi_1[imc1]-(*Phi_2)[imc2]), 2 );

      Double_t particles_distance = sqrt( dPt2 + dEta2 + dPhi2 );
      if(particles_distance < min_particles_distance){
	sel_mc_part = imc2;
	min_particles_distance = particles_distance;
	Min_dPt2_dEta2_dPhi2 = dPt2 + dEta2 + dPhi2;
      }
    }///Ends loop over MC particles

    //Flags a chosen MC as "used" and saves deltas
    if(sel_mc_part != -1){
      McFlag[sel_mc_part] = 1;
      Sum_dPt2_dEta2_dPhi2 += Min_dPt2_dEta2_dPhi2;
    }
  }///Ends loop over DATA particles
  
  
  return sqrt( Sum_dPt2_dEta2_dPhi2 );
}


void Composer(FmeSetup UserConfig){

  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Composer actived"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
  TString			TTreeName		= UserConfig.TTreeName;
  TString			IdBranch		= UserConfig.IdBranch;
  TString			PtBranch		= UserConfig.PtBranch;
  TString			EtaBranch		= UserConfig.EtaBranch;
  TString			PhiBranch		= UserConfig.PhiBranch;
  std::vector<std::string>	vMCs			= UserConfig.vMCs;
  //Double_t			ScaledPt		= UserConfig.ScaledPt;
  //Double_t			ScaledEta		= UserConfig.ScaledEta;
  //Double_t			ScaledPhi		= UserConfig.ScaledPhi;
  Int_t				GenFactor		= UserConfig.GenFactor;
  TString			DistPt			= UserConfig.DistPt;
  TString			DistEta			= UserConfig.DistEta;
  TString			DistPhi			= UserConfig.DistPhi;
  Double_t			GaussianMean		= UserConfig.GaussianMean;
  //Double_t			SDrCondition		= UserConfig.SDrCondition;
  //Double_t                    BDrCondition            = UserConfig.BDrCondition;
  //Int_t 	                MaxGenTrials            = UserConfig.MaxGenTrials;
  TString			OutGenName		= UserConfig.OutGenName;
  Float_t			MCLimit			= UserConfig.MCLimit;
  //---------------------------------------------------------------------------------
  

  //Shows with distortions will be made
  if(DistPt  == "true") std::cout<<":: Particles Pt will be smeared (gm: "<<GaussianMean<<")"<<std::endl;
  if(DistEta == "true") std::cout<<":: Particles Eta will be smeared (gm: "<<GaussianMean<<")"<<std::endl;
  if(DistPhi == "true") std::cout<<":: Particles Phi will be smeared (gm: "<<GaussianMean<<")"<<std::endl;
  
  
  ///Verifying the scale factors
  //if(ScaledPt < 0 || ScaledEta < 0 || ScaledPhi < 0){
    //std::cout<<":: ["<<ansi_yellow<<"Initials [scale_dPt/scale_dEta/scale_dPhi] -----> ["<<ScaledPt<<"/"<<ScaledEta<<"/"<<ScaledPhi<<"] @@Computing new scale factors..."<<ansi_reset<<"]"<<std::endl;
    //FindScaleFactors(vMCs, TTreeName, PtBranch, EtaBranch, PhiBranch, &ScaledPt, &ScaledEta, &ScaledPhi);
    //std::cout<<":: ["<<ansi_yellow<<"NOTE"<<ansi_reset<<Form("] Setting scale_dPt = %.3f, scale_dEta = %.3f, scale_dPhi = %.3f", ScaledPt, ScaledEta, ScaledPhi)<<std::endl;
  //}
  
  
  
  TStopwatch t2;//put outside.. does it work?!
  
  //files have to be opened here, before get number of events,
  //to avoid tree directories get overload
  TFile *inMC1 = TFile::Open((TString)vMCs[0]);
  TFile *inMC2 = TFile::Open((TString)vMCs[1]);
        
  Int_t n_sevents = get_number_of_entries((TString)vMCs[0], TTreeName);
  Int_t n_bevents = get_number_of_entries((TString)vMCs[1], TTreeName);
  Int_t max_ev_per_class = std::min(n_sevents,n_bevents);
  
  TTree *orig_sevents, *orig_bevents;
  if(MCLimit > 0){
    orig_sevents = ((TTree*)inMC1->Get(TTreeName))->CloneTree(MCLimit);
    orig_bevents = ((TTree*)inMC2->Get(TTreeName))->CloneTree(MCLimit);
  }
  else{
    orig_sevents = ((TTree*)inMC1->Get(TTreeName))->CloneTree(max_ev_per_class);
    orig_bevents = ((TTree*)inMC2->Get(TTreeName))->CloneTree(max_ev_per_class);
  }

  TList *list = new TList;
  list->Add(orig_sevents);
  list->Add(orig_bevents);
  TTree *orig_full_events = TTree::MergeTrees(list);
  list->Delete();
  inMC1->Close();
  inMC2->Close();
  
  
  Int_t oEventClass;
  std::vector<Double_t> *oParticlePt=0, *oParticleEta=0, *oParticlePhi=0;
  orig_full_events->SetBranchAddress("McFileIndex",&oEventClass);
  orig_full_events->SetBranchAddress("ParticlePt",&oParticlePt);
  orig_full_events->SetBranchAddress("ParticleEta",&oParticleEta);
  orig_full_events->SetBranchAddress("ParticlePhi",&oParticlePhi);
  Int_t nevents = orig_full_events->GetEntries();
  
  
  ///Tree to store the results from analysis
  Int_t EventClass;
  Double_t EventDist;
  std::vector<Int_t> ParticleId;
  std::vector<Double_t> ParticlePt, ParticleEta, ParticlePhi;
  TTree *sfme_tree = new TTree(TTreeName,"Events generated from FastME Generator");
  sfme_tree->SetDirectory(0);
  sfme_tree->Branch("EventClass",&EventClass,"EventClass/I");
  sfme_tree->Branch("EventDist",&EventDist,"EventDist/D");
  sfme_tree->Branch("ParticleId",&ParticleId);
  sfme_tree->Branch("ParticlePt",&ParticlePt);
  sfme_tree->Branch("ParticleEta",&ParticleEta);
  sfme_tree->Branch("ParticlePhi",&ParticlePhi);

  TTree *bfme_tree = new TTree(TTreeName,"Events generated from FastME Generator");
  bfme_tree->SetDirectory(0);
  bfme_tree->Branch("EventClass",&EventClass,"EventClass/I");
  bfme_tree->Branch("EventDist",&EventDist,"EventDist/D");
  bfme_tree->Branch("ParticleId",&ParticleId);
  bfme_tree->Branch("ParticlePt",&ParticlePt);
  bfme_tree->Branch("ParticleEta",&ParticleEta);
  bfme_tree->Branch("ParticlePhi",&ParticlePhi);

  
  //It goes over the original MCs
  std::cout<<"::   >>> Instantiating TRandom3 <<<"<<std::endl;
  TRandom3 *rg = new TRandom3();
  TH1D *gen_gaussian = new TH1D("gen_gaussian",Form("Gaussian Distribution #mu = %.2f",GaussianMean),1000,-0.05,2.05);
  
  std::cout<<"::   There are "<<nevents/2<<" events in the fountain per class"<<std::endl;
  std::cout<<"::   You resquested to multiply that by factor "<<GenFactor<<". So, generating "<<GenFactor*nevents/2<<" events per class"<<std::endl;
  std::cout<<"::   >>> Doing Data Augmentation <<<"<<std::endl;
  double distortion_factor = 1.;
  for(Int_t iev=0; iev<nevents; iev++){
    if(iev % (nevents/10) == 0)
      std::cout<<":: Remaining "<<nevents-iev<<" events"<<std::endl;

    orig_full_events->GetEntry(iev);
    Int_t nMcParticles = oParticlePt->size();

    //Each event it's distorted ntimes
    for(int itime=0; itime<GenFactor; itime++){
     
      for(Int_t imc=0; imc<nMcParticles; imc++){
	ParticleId.push_back(-1);
	
	//Pt
	if(DistPt == "true"){
	  distortion_factor = rg->Gaus(1,GaussianMean);
	  if(distortion_factor < 0){
	    while(distortion_factor < 0){
	      distortion_factor = rg->Gaus(1,GaussianMean);
	    }
	  }
          gen_gaussian->Fill(distortion_factor);
	  ParticlePt.push_back( (*oParticlePt)[imc]*distortion_factor );
	}
	else{
	  ParticlePt.push_back( (*oParticlePt)[imc] );
	}
	
	//Eta
	if(DistEta == "true"){
	  distortion_factor = rg->Gaus(1,GaussianMean);
	  if(distortion_factor < 0){
	    while(distortion_factor < 0){
	      distortion_factor = rg->Gaus(1,GaussianMean);
	    }
	  }
          gen_gaussian->Fill(distortion_factor);
	  ParticleEta.push_back( (*oParticleEta)[imc]*distortion_factor );
	}
	else{
	  ParticleEta.push_back( (*oParticleEta)[imc] );
	}

	//Phi
	if(DistPhi == "true"){
	  distortion_factor = rg->Gaus(1,GaussianMean);
	  if(distortion_factor < 0){
	    while(distortion_factor < 0){
	      distortion_factor = rg->Gaus(1,GaussianMean);
	    }
	  }
          gen_gaussian->Fill(distortion_factor);
	  ParticlePhi.push_back( (*oParticlePhi)[imc]*distortion_factor );
	}
	else{
	  ParticlePhi.push_back( (*oParticlePhi)[imc] );
	}
	
      }      
      EventClass = oEventClass;
      
      //Compute the distance to the original event
      EventDist = get_event_distance(nMcParticles, ParticlePt, ParticleEta, ParticlePhi, oParticlePt, oParticleEta, oParticlePhi);
      
      if(EventClass == 0) sfme_tree->Fill();
      else if(EventClass == 1) bfme_tree->Fill();
      else std::cout<<":: ERROR: MC type not defined!!"<<std::endl;
      ParticleId.clear();
      ParticlePt.clear();
      ParticleEta.clear();
      ParticlePhi.clear();
    }///End loop over time generation    
  }///End loop over MC events
    

  t2.Stop();
  std::cout<<"[!Event generation finished!]"<<std::endl;
  std::cout<<"[!Saving files "<<OutGenName<<"_SIG.root & "<<OutGenName<<"_BKG.root!]"<<std::endl;
  TFile *soutgen = new TFile(OutGenName + "_SIG.root","recreate");
  sfme_tree->Write();
  gen_gaussian->Write();
  soutgen->Close();

  TFile *boutgen = new TFile(OutGenName + "_BKG.root","recreate");
  bfme_tree->Write();
  gen_gaussian->Write();
  boutgen->Close();

  //cleaning memory
  delete oParticlePt;
  delete oParticleEta;
  delete oParticlePhi;
  
  return;
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}



#endif
