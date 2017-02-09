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


float gen_value(float input, float lower_limit, float upper_limit){
  if(lower_limit >= 0) return lower_limit + input*(upper_limit-lower_limit);
  else return -lower_limit + 2*input*upper_limit;
}



void Composer(FmeSetup UserConfig){

  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Composer actived"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
  TString			TTreeName		= UserConfig.TTreeName;
  TString			IdBranch		= UserConfig.IdBranch;
  TString			PtBranch		= UserConfig.PtBranch;
  TString			EtaBranch		= UserConfig.EtaBranch;
  TString			PhiBranch		= UserConfig.PhiBranch;
  std::vector<std::string>	vMCs			= UserConfig.vMCs;
  TString			SetFlavorConstraint	= UserConfig.SetFlavorConstraint;
  Double_t			ScaledPt		= UserConfig.ScaledPt;
  Double_t			ScaledEta		= UserConfig.ScaledEta;
  Double_t			ScaledPhi		= UserConfig.ScaledPhi;
  Int_t				GenNEv			= UserConfig.GenNEv;
  Double_t			SDrCondition		= UserConfig.SDrCondition;
  Double_t                      BDrCondition            = UserConfig.BDrCondition;
  Int_t 	                MaxGenTrials            = UserConfig.MaxGenTrials;
  TString			OutGenName		= UserConfig.OutGenName;
  Float_t			MCLimit			= UserConfig.MCLimit;
  //---------------------------------------------------------------------------------
  
  
  
  
  ///Verifying the scale factors
  if(ScaledPt < 0 || ScaledEta < 0 || ScaledPhi < 0){
    std::cout<<":: ["<<ansi_yellow<<"Initials [scale_dPt/scale_dEta/scale_dPhi] -----> ["<<ScaledPt<<"/"<<ScaledEta<<"/"<<ScaledPhi<<"] @@Computing new scale factors..."<<ansi_reset<<"]"<<std::endl;
    FindScaleFactors(vMCs, TTreeName, PtBranch, EtaBranch, PhiBranch, &ScaledPt, &ScaledEta, &ScaledPhi);
    std::cout<<":: ["<<ansi_yellow<<"NOTE"<<ansi_reset<<Form("] Setting scale_dPt = %.3f, scale_dEta = %.3f, scale_dPhi = %.3f", ScaledPt, ScaledEta, ScaledPhi)<<std::endl;
  }
  
  
  
  TStopwatch t2;//put outside.. does it work?!
        
  ///Addresses the MC branches to be used
  TFile *inMC1 = TFile::Open((TString)vMCs[0]);
  TTreeReader tread1(TTreeName,inMC1);  
  TTreeReaderValue<int>    	McType1(tread1, "McFileIndex");
  TTreeReaderArray<double>	McPt1(tread1, PtBranch);
  TTreeReaderArray<double>	McEta1(tread1, EtaBranch);
  TTreeReaderArray<double>	McPhi1(tread1, PhiBranch);

  ///Addresses the MC branches to be used
  TFile *inMC2 = TFile::Open((TString)vMCs[1]);
  TTreeReader tread2(TTreeName,inMC2);  
  TTreeReaderValue<int>    	McType2(tread2, "McFileIndex");
  TTreeReaderArray<double>	McPt2(tread2, PtBranch);
  TTreeReaderArray<double>	McEta2(tread2, EtaBranch);
  TTreeReaderArray<double>	McPhi2(tread2, PhiBranch);

  
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

  
  std::cout<<":: Defining parameters for generation..."<<std::endl;
  Int_t nMonteCarlo = tread1.GetEntries(true);
  tread1.SetEntry(0);
  Int_t nMcParticles = McPt1.GetSize();
  Int_t nDataParticles = nMcParticles;
  //std::vector<float> lower_limit, upper_limit;
  //for(int b=0; b<nMcParticles*3; b++){
    //lower_limit.push_back(99);
    //upper_limit.push_back(-99);
  //}
  TH1 *for_pt  = new TH1D("for_pt","for_pt",200,0,1000);
  TH1 *for_eta = new TH1D("for_eta","for_eta",20,-10,10);
  float phi_low_limit = 99, phi_high_limit = -99;
  for(Int_t mc = 0; mc < nMonteCarlo; ++mc){
    tread1.SetEntry(mc); ///Move on MC loop	

    for(Int_t imc = 0; imc < nMcParticles; ++imc){
      //if( McPt1[imc]  < lower_limit[3*imc] ) lower_limit[3*imc] = McPt1[imc];
      //if( McEta1[imc] < lower_limit[3*imc+1] ) lower_limit[3*imc+1] = McEta1[imc];
      //if( McPhi1[imc] < lower_limit[3*imc+2] ) lower_limit[3*imc+2] = McPhi1[imc];
      //if( McPt1[imc]  > upper_limit[3*imc] ) upper_limit[3*imc] = McPt1[imc];
      //if( McEta1[imc] > upper_limit[3*imc+1] ) upper_limit[3*imc+1] = McEta1[imc];
      //if( McPhi1[imc] > upper_limit[3*imc+2] ) upper_limit[3*imc+2] = McPhi1[imc];
      for_pt->Fill( McPt1[imc] );
      for_eta->Fill( McEta1[imc] );
      phi_low_limit  = (McPhi1[imc] < phi_low_limit)? McPhi1[imc] : phi_low_limit;
      phi_high_limit = (McPhi1[imc] > phi_high_limit)? McPhi1[imc] : phi_high_limit;
    }
  }
  std::cout<<":: ---------- Extracted Parameters -----------"<<std::endl;
  //for(int i=0; i<nMcParticles; i++){
    //std::cout<<":: Particle "<<i<<std::endl;
    //std::cout<<":: Lower Pt: "<<lower_limit[3*i]<<"\t\tUpper Pt: "<<upper_limit[3*i]<<std::endl;
    //std::cout<<":: Lower Eta: "<<lower_limit[3*i+1]<<"\t\tUpper Eta: "<<upper_limit[3*i+1]<<std::endl;
    //std::cout<<":: Lower Phi: "<<lower_limit[3*i+2]<<"\t\tUpper Phi: "<<upper_limit[3*i+2]<<std::endl;
  //}
  float pt_mean  = for_pt->GetMean();
  float pt_rms = for_pt->GetRMS();
  float eta_mean  = for_eta->GetMean();
  float eta_rms = for_eta->GetRMS();
  std::cout<<":: PtMean:   "<<pt_mean<<"\tPtRMS:   "<<pt_rms<<std::endl;
  std::cout<<":: EtaMean:  "<<eta_mean<<"\tEtaRMS:  "<<eta_rms<<std::endl;
  std::cout<<":: LowerPhi: "<<phi_low_limit<<"\tUpperPhi: "<<phi_high_limit<<std::endl;
  std::cout<<":: -----------------------------"<<std::endl;  
  
  int SAcpdEvents = 0, BAcpdEvents = 0, itrial = 0;
  std::vector<TRandom3> GenEngine(nMcParticles*3);
  for(int g=0; g<nMcParticles*3; g++) GenEngine[g].SetSeed(g*100);
  std::vector<double> DataPt, DataEta, DataPhi;
  for(int b=0; b<nMcParticles; b++){
    DataPt.push_back(0);
    DataEta.push_back(0);
    DataPhi.push_back(0);
  }

  std::cout<<":: "<<ansi_yellow<<">>>>> SAMPLING <<<<<"<<ansi_reset<<std::endl;
  while(SAcpdEvents < GenNEv || BAcpdEvents < GenNEv){
    if(itrial > MaxGenTrials) break;
    itrial++;
  
    //std::cout<<":: Event gen..."<<std::endl;
    for(Int_t idt=0; idt<nMcParticles; idt++){
      DataPt[idt]  = fabs( GenEngine[3*idt].Gaus(pt_mean,pt_rms) );
      DataEta[idt] = GenEngine[3*idt+1].Gaus(eta_mean,eta_rms);
      DataPhi[idt] = gen_value( GenEngine[3*idt+2].Rndm(), phi_low_limit, phi_high_limit );
      //DataPt[idt] = gen_value(GenEngine[3*idt].Rndm(), lower_limit[3*idt], upper_limit[3*idt]);
      //DataEta[idt] = gen_value(GenEngine[3*idt+1].Rndm(), lower_limit[3*idt+1], upper_limit[3*idt+1]);
      //DataPhi[idt] = gen_value(GenEngine[3*idt+2].Rndm(), lower_limit[3*idt+2], upper_limit[3*idt+2]);

      //std::cout<<":: Pt: "<<DataPt[idt]<<"\tEta: "<<DataEta[idt]<<"\tPhi: "<<DataPhi[idt]<<std::endl;
    }
    
    float min_distance_gen1 = 1.e15, min_distance_gen2 = 1.e15;
    for(Int_t mc = 0; mc < (Int_t)tread1.GetEntries(true); ++mc){
      if(mc > (Int_t)MCLimit) break;
      tread1.SetEntry(mc); ///Move on MC loop
      std::vector<int> McFlag(nMcParticles,0);
	  
      //Loop over the particles in the data event
      Double_t Min_dPt2_dEta2_dPhi2 = 0, Sum_dPt2_dEta2_dPhi2 = 0;
      for(int idt = 0; idt < nDataParticles; ++idt){
	Double_t min_particles_distance = 1.e15;
	Int_t sel_mc_part = -1;
	    
	//Loop over the particles in the MC event
	for(Int_t imc = 0; imc < nMcParticles; ++imc){
	  ///Avoid MC flaged as "used"
	  if(McFlag[imc] == 1) continue;
	      
	  ///Compute preliminary particles distance
	  Double_t dPt2  = pow( (DataPt[idt]-McPt1[imc])*ScaledPt, 2 );	      
	  Double_t dEta2 = pow( (DataEta[idt]-McEta1[imc])*ScaledEta, 2 );
	  Double_t dPhi2 = pow( (DataPhi[idt]-McPhi1[imc])*ScaledPhi, 2 );

	  Double_t particles_distance = sqrt( dPt2 + dEta2 + dPhi2 );
	  //if( Verbose == 3 ) std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<sMcId[imc]<<"   part_dist: "<<particles_distance<<std::endl;
	  if(particles_distance < min_particles_distance){
	    sel_mc_part = imc;
	    min_particles_distance = particles_distance;
	    Min_dPt2_dEta2_dPhi2 = dPt2 + dEta2 + dPhi2;
	  }
	}///Ends loop over gen particles

	//Flags a chosen MC as "used" and saves deltas
	if(sel_mc_part != -1){
	  McFlag[sel_mc_part] = 0;
	  //Start to sum for final event distance
	  Sum_dPt2_dEta2_dPhi2 += Min_dPt2_dEta2_dPhi2;
	}
      }///Ends loop over DATA particles
	
      //Computes the final event distance
      if(sqrt(Sum_dPt2_dEta2_dPhi2) < min_distance_gen1){
	min_distance_gen1 = sqrt(Sum_dPt2_dEta2_dPhi2);
      }
    }///End loop over MC events

    for(Int_t mc = 0; mc < (Int_t)tread2.GetEntries(true); ++mc){
      if(mc > (Int_t)MCLimit) break;
      tread2.SetEntry(mc); ///Move on MC loop
      std::vector<int> McFlag(nMcParticles,0);
	  
      //Loop over the particles in the data event
      Double_t Min_dPt2_dEta2_dPhi2 = 0, Sum_dPt2_dEta2_dPhi2 = 0;
      for(int idt = 0; idt < nDataParticles; ++idt){
	Double_t min_particles_distance = 1.e15;
	Int_t sel_mc_part = -1;
	    
	//Loop over the particles in the MC event
	for(Int_t imc = 0; imc < nMcParticles; ++imc){
	  ///Avoid MC flaged as "used"
	  if(McFlag[imc] == 1) continue;
	      
	  ///Compute preliminary particles distance
	  Double_t dPt2  = pow( (DataPt[idt]-McPt2[imc])*ScaledPt, 2 );	      
	  Double_t dEta2 = pow( (DataEta[idt]-McEta2[imc])*ScaledEta, 2 );
	  Double_t dPhi2 = pow( (DataPhi[idt]-McPhi2[imc])*ScaledPhi, 2 );

	  Double_t particles_distance = sqrt( dPt2 + dEta2 + dPhi2 );
	  //if( Verbose == 3 ) std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<sMcId[imc]<<"   part_dist: "<<particles_distance<<std::endl;
	  if(particles_distance < min_particles_distance){
	    sel_mc_part = imc;
	    min_particles_distance = particles_distance;
	    Min_dPt2_dEta2_dPhi2 = dPt2 + dEta2 + dPhi2;
	  }
	}///Ends loop over gen particles

	//Flags a chosen MC as "used" and saves deltas
	if(sel_mc_part != -1){
	  McFlag[sel_mc_part] = 0;
	  //Start to sum for final event distance
	  Sum_dPt2_dEta2_dPhi2 += Min_dPt2_dEta2_dPhi2;
	}
      }///Ends loop over DATA particles
	
      //Computes the final event distance
      if(sqrt(Sum_dPt2_dEta2_dPhi2) < min_distance_gen2){
	min_distance_gen2 = sqrt(Sum_dPt2_dEta2_dPhi2);
      }
    }///End loop over MC2 events
    
    float disc = min_distance_gen2/(min_distance_gen1 + min_distance_gen2);
    EventClass = (min_distance_gen1 < min_distance_gen2)? *McType1 : *McType2;
    if( disc < BDrCondition || disc > SDrCondition){
      EventDist = disc;
      for(Int_t idt=0; idt<nMcParticles; idt++){
	ParticleId.push_back(-1);
	ParticlePt.push_back(DataPt[idt]);
	ParticleEta.push_back(DataEta[idt]);
	ParticlePhi.push_back(DataPhi[idt]);
      }

      if(SAcpdEvents <= GenNEv && EventClass == *McType1){ sfme_tree->Fill(); SAcpdEvents++; }
      if(BAcpdEvents <= GenNEv && EventClass == *McType2){ bfme_tree->Fill(); BAcpdEvents++; }
      ParticleId.clear();
      ParticlePt.clear();
      ParticleEta.clear();
      ParticlePhi.clear();
    }
      
    if(itrial % 50 == 0) 
	std::cout<<":: Found S: "<<ansi_cyan<<SAcpdEvents<<ansi_reset<<" and B: "<<ansi_red<<BAcpdEvents<<ansi_reset<<" events in "<<itrial<<" trials..."<<std::endl;
  }///End loop over trials
    

  t2.Stop();
  std::cout<<"[!Event generation finished!]"<<std::endl;
  std::cout<<"[!Saving files "<<OutGenName<<"_sig.root & "<<OutGenName<<"_bkg.root!]"<<std::endl;
  TFile *soutgen = new TFile(OutGenName + "_SIG.root","recreate");
  sfme_tree->Write();
  soutgen->Close();

  TFile *boutgen = new TFile(OutGenName + "_BKG.root","recreate");
  bfme_tree->Write();
  boutgen->Close();

  
  return;
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}



#endif
