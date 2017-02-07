///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::[ Cartographer  - Compute Events Distance ]::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Generator_h
#define Generator_h


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


/*
///Compute scale factors to normalize the deltas
void FindScaleFactors(std::vector<std::string> vMCs, TString TTreeName, TString PtBranch, TString EtaBranch, 
		      TString PhiBranch, Double_t *ScaledPt, Double_t *ScaledEta, Double_t *ScaledPhi){
  
  gROOT->SetBatch();
  TCanvas *temp = new TCanvas();
  Double_t pt_sum = 0, eta_sum = 0, phi_sum = 0;
  Int_t total = vMCs.size();
  for(Int_t isample = 0; isample < total; isample++){

    TFile *ftmp = TFile::Open( (TString)vMCs[isample] );
    TTree *ttmp = (TTree*)ftmp->Get( TTreeName );
    
    if( *ScaledPt < 0){
      TString draw_pt = PtBranch + " >> stackpt";
      ttmp->Draw(draw_pt);
      TH1D *stackpt = (TH1D*)gDirectory->Get("stackpt");
      if( *ScaledPt == -1 ) pt_sum += stackpt->GetMean();
      if( *ScaledPt == -2 ) pt_sum += stackpt->GetBinCenter( stackpt->GetMaximumBin() );
    }
    if( *ScaledEta < 0){
      TString draw_eta = EtaBranch+" >> stacketa";
      ttmp->Draw(draw_eta);
      TH1D *stacketa = (TH1D*)gDirectory->Get("stacketa");
      if( *ScaledEta == -1 ) eta_sum += fabs( stacketa->GetMean() );//Should be used only in case you are in region shifted from 0!
      if( *ScaledEta == -2 ) eta_sum += fabs( stacketa->GetBinCenter(stacketa->GetMinimum()) );
    }
    if( *ScaledPhi < 0){
      TString draw_phi = PhiBranch+" >> stackphi";
      ttmp->Draw(draw_phi);
      TH1D *stackphi = (TH1D*)gDirectory->Get("stackphi");
      if( *ScaledPhi == -1 ) phi_sum += fabs( stackphi->GetMean() );//Should be used only in case you are in region shifted from 0!
      if( *ScaledPhi == -2 ) phi_sum += fabs( stackphi->GetBinCenter(stackphi->GetMinimum()) );
    }
   
    ftmp->Close();
  }
  
  if( *ScaledPt  < 0 )  *ScaledPt  = total/pt_sum;
  if( *ScaledEta < 0 )  *ScaledEta = total/eta_sum;
  if( *ScaledPhi < 0 )  *ScaledPhi = total/phi_sum;
  
  temp->Close();
  
  return;
}
*/




///######################################## FastME Main Function: Finds the closest events and store the distances ###################################################
void Generator(FmeSetup UserConfig){

  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Getting User Configuration"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
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
  Double_t			DrCondition		= UserConfig.DrCondition;
  TString			OutGenName		= UserConfig.OutGenName;
  Int_t				Verbose			= UserConfig.Verbose;
  //---------------------------------------------------------------------------------
  
  
  
  
  ///Verifying the scale factors
  if(ScaledPt < 0 || ScaledEta < 0 || ScaledPhi < 0){
    std::cout<<":: ["<<ansi_yellow<<"Initials [scale_dPt/scale_dEta/scale_dPhi] -----> ["<<ScaledPt<<"/"<<ScaledEta<<"/"<<ScaledPhi<<"] @@Computing new scale factors..."<<ansi_reset<<"]"<<std::endl;
    FindScaleFactors(vMCs, TTreeName, PtBranch, EtaBranch, PhiBranch, &ScaledPt, &ScaledEta, &ScaledPhi);
    std::cout<<":: ["<<ansi_yellow<<"NOTE"<<ansi_reset<<Form("] Setting scale_dPt = %.3f, scale_dEta = %.3f, scale_dPhi = %.3f", ScaledPt, ScaledEta, ScaledPhi)<<std::endl;
  }
  
  
  
  TStopwatch t2;//put outside.. does it work?!
        
  ///Addresses the MC branches to be used
  TFile *inMC = TFile::Open((TString)vMCs[0]);
  TTreeReader tread(TTreeName,inMC);  
  TTreeReaderValue<int>    	McType(tread, "McFileIndex");
  TTreeReaderArray<double>	McPt(tread, PtBranch);
  TTreeReaderArray<double>	McEta(tread, EtaBranch);
  TTreeReaderArray<double>	McPhi(tread, PhiBranch);


  ///Tree to store the results from analysis
  Int_t EventClass;
  Double_t EventDist;
  std::vector<Int_t> ParticleId;
  std::vector<Double_t> ParticlePt, ParticleEta, ParticlePhi;
  TTree *fme_tree = new TTree(TTreeName,"Events generated from FastME Generator");
  fme_tree->SetDirectory(0);
  fme_tree->Branch("EventClass",&EventClass,"EventClass/I");
  fme_tree->Branch("EventDist",&EventDist,"EventDist/D");
  fme_tree->Branch("ParticleId",&ParticleId);
  fme_tree->Branch("ParticlePt",&ParticlePt);
  fme_tree->Branch("ParticleEta",&ParticleEta);
  fme_tree->Branch("ParticlePhi",&ParticlePhi);

  
  std::cout<<":: Defining the limits for generation..."<<std::endl;
  Int_t nMonteCarlo = tread.GetEntries(true);
  tread.SetEntry(0);
  Int_t nMcParticles = McPt.GetSize();
  Int_t nDataParticles = nMcParticles;
  std::vector<float> lower_limit, upper_limit;
  for(int b=0; b<nMcParticles*3; b++){
    lower_limit.push_back(99);
    upper_limit.push_back(-99);
  }
  for(Int_t mc = 0; mc < nMonteCarlo; ++mc){
    tread.SetEntry(mc); ///Move on MC loop	

    for(Int_t imc = 0; imc < nMcParticles; ++imc){
      if( McPt[imc]  < lower_limit[3*imc] ) lower_limit[3*imc] = McPt[imc];
      if( McEta[imc] < lower_limit[3*imc+1] ) lower_limit[3*imc+1] = McEta[imc];
      if( McPhi[imc] < lower_limit[3*imc+2] ) lower_limit[3*imc+2] = McPhi[imc];
      if( McPt[imc]  > upper_limit[3*imc] ) upper_limit[3*imc] = McPt[imc];
      if( McEta[imc] > upper_limit[3*imc+1] ) upper_limit[3*imc+1] = McEta[imc];
      if( McPhi[imc] > upper_limit[3*imc+2] ) upper_limit[3*imc+2] = McPhi[imc];
    }
  }
  std::cout<<":: ---------- Limits -----------"<<std::endl;
  for(int i=0; i<nMcParticles; i++){
    std::cout<<":: Particle "<<i<<std::endl;
    std::cout<<":: Lower Pt: "<<lower_limit[3*i]<<"\t\tUpper Pt: "<<upper_limit[3*i]<<std::endl;
    std::cout<<":: Lower Eta: "<<lower_limit[3*i+1]<<"\t\tUpper Eta: "<<upper_limit[3*i+1]<<std::endl;
    std::cout<<":: Lower Phi: "<<lower_limit[3*i+2]<<"\t\tUpper Phi: "<<upper_limit[3*i+2]<<std::endl;
  }
  std::cout<<":: -----------------------------"<<std::endl;  
  
  int AcpdEvents = 0, itrial = 0;
  std::vector<TRandom3> GenEngine(nMcParticles*3);
  for(int g=0; g<nMcParticles*3; g++) GenEngine[g].SetSeed(g*100);
  std::vector<double> DataPt, DataEta, DataPhi;
  for(int b=0; b<nMcParticles; b++){
    DataPt.push_back(0);
    DataEta.push_back(0);
    DataPhi.push_back(0);
  }

  float min_distance_gen = 1.e15;
  std::cout<<":: Starting sampling..."<<std::endl;
  while(AcpdEvents < GenNEv){
    itrial++;
    float average_dr = 0;
  
    //std::cout<<":: Event gen..."<<std::endl;
    for(Int_t idt=0; idt<nMcParticles; idt++){
      DataPt[idt] = gen_value(GenEngine[3*idt].Rndm(), lower_limit[3*idt], upper_limit[3*idt]);
      DataEta[idt] = gen_value(GenEngine[3*idt+1].Rndm(), lower_limit[3*idt+1], upper_limit[3*idt+1]);
      DataPhi[idt] = gen_value(GenEngine[3*idt+2].Rndm(), lower_limit[3*idt+2], upper_limit[3*idt+2]);

      //std::cout<<":: Pt: "<<DataPt[idt]<<"\tEta: "<<DataEta[idt]<<"\tPhi: "<<DataPhi[idt]<<std::endl;
    }
    
    for(Int_t mc = 0; mc < nMonteCarlo; ++mc){
      tread.SetEntry(mc); ///Move on MC loop
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
	  Double_t dPt2  = pow( (DataPt[idt]-McPt[imc])*ScaledPt, 2 );	      
	  Double_t dEta2 = pow( (DataEta[idt]-McEta[imc])*ScaledEta, 2 );
	  Double_t dPhi2 = pow( (DataPhi[idt]-McPhi[imc])*ScaledPhi, 2 );

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
	
/*	  
      //Computes the final event distance
      if(sqrt(Sum_dPt2_dEta2_dPhi2) < DrCondition){
	AcpdEvents++;
	EventClass = *McType;
        EventDist = sqrt(Sum_dPt2_dEta2_dPhi2);
	for(Int_t idt=0; idt<nMcParticles; idt++){
	  ParticleId.push_back(-1);
	  ParticlePt.push_back(DataPt[idt]);
	  ParticleEta.push_back(DataEta[idt]);
	  ParticlePhi.push_back(DataPhi[idt]);
	}
	fme_tree->Fill();
	ParticleId.clear();
	ParticlePt.clear();
	ParticleEta.clear();
	ParticlePhi.clear();
      }
*/    
      //if(sqrt(Sum_dPt2_dEta2_dPhi2) < min_distance_gen) min_distance_gen = sqrt(Sum_dPt2_dEta2_dPhi2);
      average_dr += sqrt(Sum_dPt2_dEta2_dPhi2)/float(nMonteCarlo);
    }///End loop over MC evens
      if(fabs(average_dr - DrCondition) < 3){
        AcpdEvents++;
        EventClass = *McType;
        EventDist = average_dr;
        for(Int_t idt=0; idt<nMcParticles; idt++){
          ParticleId.push_back(-1);
          ParticlePt.push_back(DataPt[idt]);
          ParticleEta.push_back(DataEta[idt]);
          ParticlePhi.push_back(DataPhi[idt]);
        }
        fme_tree->Fill();
        ParticleId.clear();
        ParticlePt.clear();
        ParticleEta.clear();
        ParticlePhi.clear();

      }
      if(itrial % 100 == 0) std::cout<<":: Found "<<AcpdEvents<<" events in "<<itrial<<" trials -- Min.Dist so far "<<min_distance_gen<<std::endl;
  }///End loop over trials
      

  t2.Stop();
  std::cout<<"[!Event generation finished!]"<<std::endl;
  std::cout<<"[!Saving file "<<OutGenName<<"!]"<<std::endl;
  TFile *outgen = new TFile(OutGenName,"recreate");
  fme_tree->Write();
  outgen->Close();
  
  return;
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}



#endif
