///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::[ Cartographer  - Compute Events Distance ]::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Cartographer_h
#define Cartographer_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/Scaler.h"

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

///Headers to TProcPool
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TProcPool.h>
#include <PoolUtils.h>
#include <TSystem.h>




///######################################## FastME Main Function: Finds the closest events and store the distances ###################################################
TTree *Cartographer(FmeSetup UserConfig){

  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Getting User Configuration"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
  std::vector<std::string>	vDatas			= UserConfig.vDatas;
  Int_t				DTLimit			= UserConfig.DTLimit;
  TString			TTreeName		= UserConfig.TTreeName;
  TString			IdBranch		= UserConfig.IdBranch;
  TString			PtBranch		= UserConfig.PtBranch;
  TString			EtaBranch		= UserConfig.EtaBranch;
  TString			PhiBranch		= UserConfig.PhiBranch;
  std::vector<std::string>	vMCs			= UserConfig.vMCs;
  UInt_t			NCores			= UserConfig.NCores;
  TString			PhSDrMethod		= UserConfig.PhSDrMethod;
  TString			SetFlavorConstraint	= UserConfig.SetFlavorConstraint;
  Float_t			MCLimit			= UserConfig.MCLimit;
  Double_t			ScaledPt		= UserConfig.ScaledPt;
  Double_t			ScaledEta		= UserConfig.ScaledEta;
  Double_t			ScaledPhi		= UserConfig.ScaledPhi;
  Int_t				Verbose			= UserConfig.Verbose;
  //---------------------------------------------------------------------------------
  
  
  
  
  ///Verifying the scale factors
  if(ScaledPt < 0 || ScaledEta < 0 || ScaledPhi < 0){
    std::cout<<":: ["<<ansi_yellow<<"Initials [scale_dPt/scale_dEta/scale_dPhi] -----> ["<<ScaledPt<<"/"<<ScaledEta<<"/"<<ScaledPhi<<"] @@Computing new scale factors..."<<ansi_reset<<"]"<<std::endl;
    FindScaleFactors(vMCs, TTreeName, PtBranch, EtaBranch, PhiBranch, &ScaledPt, &ScaledEta, &ScaledPhi);
    std::cout<<":: ["<<ansi_yellow<<"NOTE"<<ansi_reset<<Form("] Setting scale_dPt = %.3f, scale_dEta = %.3f, scale_dPhi = %.3f", ScaledPt, ScaledEta, ScaledPhi)<<std::endl;
  }
  
  
  
  ///Mappers Definition
  //===================================================================================================================================================================
  
  ///Definition of Minimum Particles Distance Method to map the events
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  auto MPD = [IdBranch, PtBranch, EtaBranch, PhiBranch, vDatas, TTreeName, DTLimit, MCLimit, ScaledPt, ScaledEta,
	      ScaledPhi, SetFlavorConstraint, Verbose](TTreeReader &tread)->TObject* {
		     
    std::cout<<"Starting MPD..."<<std::endl;
    TStopwatch t2;//put outside.. does it work?!
        
    ///Addresses the MC branches to be used
    TTreeReaderValue<int>    	McType(tread, "McFileIndex");
    TTreeReaderArray<int>	McId(tread, IdBranch);
    TTreeReaderArray<double>	McPt(tread, PtBranch);
    TTreeReaderArray<double>	McEta(tread, EtaBranch);
    TTreeReaderArray<double>	McPhi(tread, PhiBranch);


    ///Tree to store the results from analysis
    Int_t DataFile;
    std::vector<double> MinDistance;
    std::vector<int> PairedMc, McFile;
    TTree *fme_tree = new TTree("fme_tree","From FastME Phase Space Analysis");
    fme_tree->SetDirectory(0);
    fme_tree->Branch("DataFile",&DataFile,"DataFile/I");
    fme_tree->Branch("MinDistance","std::vector<double>",&MinDistance);
    fme_tree->Branch("PairedMc","std::vector<int>",&PairedMc);
    fme_tree->Branch("McFile","std::vector<int>",&McFile);


    ///Loop over the different data files
    //std::cout<<"Over the loop now..."<<std::endl;
    Int_t nDataFiles = vDatas.size();
    for(Int_t i_data_file = 0; i_data_file < nDataFiles; i_data_file++){
      DataFile = i_data_file;
      MinDistance.clear();
      PairedMc.clear();
      McFile.clear();


      TFile *fData = TFile::Open( (TString)vDatas[i_data_file] );
      TTreeReader refReader(TTreeName, fData);
      TTreeReaderArray<int>      DataId(refReader, IdBranch);
      TTreeReaderArray<double>   DataPt(refReader, PtBranch);
      TTreeReaderArray<double>   DataEta(refReader, EtaBranch);
      TTreeReaderArray<double>   DataPhi(refReader, PhiBranch);

    
      ///Loop on Data events
      Int_t nDataEvents = refReader.GetEntries(true);
      if(DTLimit >= 1 && DTLimit < nDataEvents) nDataEvents = DTLimit;
      MinDistance.reserve(nDataEvents);
      PairedMc.reserve(nDataEvents);
      McFile.reserve(nDataEvents);
      //std::cout<<"NData: "<<nDataEvents<<std::endl;
      for(Int_t dt = 0; dt < nDataEvents; dt++){
	//std::cout<<"Pulling data entry: "<<dt<<std::endl;
	refReader.SetEntry(dt); ///Move on Data loop     
	Int_t nDataParticles = DataPt.GetSize();
	//std::cout<<"Pulled successfully..."<<std::endl;                                                         
      
	if( Verbose != 0 && ((dt!= 0 && nDataEvents > 10 && dt%(nDataEvents/10) == 0) || (nDataEvents-dt) == 1) ){
	  std::string infos = Form("%i/%i/%i", i_data_file, *McType, nDataEvents-dt);
	  std::cout<<"\r:: DataFile/MC/Remaining: "<<infos<<std::endl;
	}


	Double_t min_distance_Min = 1.e15;
	Int_t imc_min = -1;
	Int_t f_type=-99;
	Int_t nMonteCarlo = tread.GetEntries(true);
	if(MCLimit >= 1 && MCLimit <= nMonteCarlo) nMonteCarlo = MCLimit;
	if(MCLimit > 0 && MCLimit < 1 && (Int_t)(MCLimit*nMonteCarlo) != 0) nMonteCarlo = MCLimit*nMonteCarlo;	
        
	//std::cout<<"Going over MC loop..."<<std::endl;
	for(Int_t mc = 0; mc < nMonteCarlo; ++mc){
	  tread.SetEntry(mc); ///Move on MC loop
	  //std::cout<<"Successfully pulled MC..."<<std::endl;
	  ///Avoid different final state comparison
	  Int_t nMcParticles = McPt.GetSize();
	  if( nDataParticles != nMcParticles ) continue;
	  std::vector<int> McFlag(nMcParticles,0);
	  
	  //Loop over the particles in the data event
	  //std::cout<<"Going over data objects..."<<std::endl;
	  Int_t n_matched_particles = 0;
	  Double_t Min_dPt2_dEta2_dPhi2 = 0, Sum_dPt2_dEta2_dPhi2 = 0;
	  for(int idt = 0; idt < nDataParticles; ++idt){
	    Double_t min_particles_distance = 1.e15;
	    Int_t sel_mc_part = -1;
	    
	    //Loop over the particles in the MC event
	    //std::cout<<"VectorSize: "<<sMcId.size()<<std::endl;
	    for(Int_t imc = 0; imc < nMcParticles; ++imc){
	      //std::cout<<"Data "<<idt<<"--- Mc "<<imc<<std::endl;

	      ///Avoid MC flaged as "used"
	      if(McFlag[imc] == 1) continue;
	      
	      ///Avoid different Data-MC particles comparison in case of flavor constraint
	      if(SetFlavorConstraint == "true" && DataId[idt] != McId[imc]) continue;

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
	    }///Ends loop over MC particles

	    //Flags a chosen MC as "used" and saves deltas
	    //std::cout<<"Selected MC "<<sel_mc_part<<std::endl;
	    if(sel_mc_part != -1){
	      McFlag[sel_mc_part] = 1;
	      ++n_matched_particles;
	      //Start to sum for final event distance
	      Sum_dPt2_dEta2_dPhi2 += Min_dPt2_dEta2_dPhi2;
	    }
	  }///Ends loop over DATA particles
	
	  
	  //Computes the final event distance
	  if(sqrt(Sum_dPt2_dEta2_dPhi2) < min_distance_Min && n_matched_particles == nDataParticles){
	    min_distance_Min = sqrt(Sum_dPt2_dEta2_dPhi2);
	    imc_min = mc;
	    f_type = *McType;	
	  }
	  //if( Verbose > 2 ) std::cout<<"Event_distance (MinDr) = "<<event_distance_Min<<std::endl;  
        }///End loop over MC evens

      
	///Stores the minimum distances found for each data event
	MinDistance.push_back(min_distance_Min);
	McFile.push_back(f_type);
	PairedMc.push_back(imc_min);
	//if( Verbose > 1 ) std::cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Min<<std::endl;
      }///End loop over Data events
  
      fme_tree->Fill();
      delete fData;
    }///Ends the loop over data files
    

    t2.Stop();
    std::cout<<"Finished MPD..."<<std::endl;
    return fme_tree;
  };
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
//										|
//										|
//						------------------------------- + -------------------------------
//										|
//										|
										
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ///Definition of Mean of Similar Combinations Method to map the events
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  auto MSC = [IdBranch, PtBranch, EtaBranch, PhiBranch, vDatas, TTreeName, DTLimit, MCLimit, ScaledPt, ScaledEta, ScaledPhi, SetFlavorConstraint, Verbose](TTreeReader &tread2)->TObject* {
		     
    std::cout<<"Starting MSC..."<<std::endl;
    TStopwatch t2;//put outside.. does it work?!
        
    ///Addresses the MC branches to be used
    TTreeReaderValue<int>    	McType(tread2, "McFileIndex");
    TTreeReaderArray<int>	McId(tread2, IdBranch);
    TTreeReaderArray<double>	McPt(tread2, PtBranch);
    TTreeReaderArray<double>	McEta(tread2, EtaBranch);
    TTreeReaderArray<double>	McPhi(tread2, PhiBranch);


    ///Tree to store the results from analysis
    Int_t DataFile;
    std::vector<double> MinDistance;
    std::vector<int> PairedMc, McFile;
    TTree *fme_tree = new TTree("fme_tree","From FastME Phase Space Analysis");
    fme_tree->SetDirectory(0);
    fme_tree->Branch("DataFile",&DataFile,"DataFile/I");
    fme_tree->Branch("MinDistance","std::vector<double>",&MinDistance);
    fme_tree->Branch("PairedMc","std::vector<int>",&PairedMc);
    fme_tree->Branch("McFile","std::vector<int>",&McFile);
 

    ///Loop over the different data files
    for(Int_t i_data_file = 0; i_data_file < (Int_t)vDatas.size(); i_data_file++){
      DataFile = i_data_file;
      MinDistance.clear();
      PairedMc.clear();
      McFile.clear();

      TFile *fData = TFile::Open( (TString)vDatas[i_data_file] );
      TTreeReader refReader(TTreeName, fData);
      TTreeReaderArray<int>      DataId(refReader, IdBranch);
      TTreeReaderArray<double>   DataPt(refReader, PtBranch);
      TTreeReaderArray<double>   DataEta(refReader, EtaBranch);
      TTreeReaderArray<double>   DataPhi(refReader, PhiBranch);

    
      ///Loop on Data events
      Int_t nDataEvents = refReader.GetEntries(true);
      if(DTLimit >= 1 && DTLimit < nDataEvents) nDataEvents = DTLimit;
      //std::cout<<"Going over data file "<<i_data_file<<" with "<<nDataEvents<<" events"<<std::endl;
      for(Int_t dt = 0; dt < nDataEvents; ++dt){
	refReader.SetEntry(dt); ///Move on Data loop
      
	if( Verbose != 0 && ((dt!= 0 && nDataEvents > 10 && dt%(nDataEvents/10) == 0) || (nDataEvents-dt) == 1) ){
	  std::string infos = Form("%i/%i/%i", i_data_file, *McType, nDataEvents-dt);
	  std::cout<<":: DataFile/MC/Remaining: "<<infos<<std::endl;
	}


	Double_t min_distance_Med = 1.e15;
	Int_t imc_min = -1;
	Int_t f_type = -99;
	Int_t nMonteCarlo = tread2.GetEntries(true);
	if(MCLimit >= 1 && MCLimit <= nMonteCarlo) nMonteCarlo = MCLimit;
	if(MCLimit > 0 && MCLimit < 1 && (Int_t)(MCLimit*nMonteCarlo) != 0) nMonteCarlo = MCLimit*nMonteCarlo;	
        
        //std::cout<<"Throwing event "<<dt<<" into "<<nMonteCarlo<<" MC events"<<std::endl;
	for(Int_t mc = 0; mc < nMonteCarlo; ++mc){
	  tread2.SetEntry(mc); ///Move on MC loop
          Int_t nDataParticles = DataPt.GetSize();
	  Int_t nMcParticles   = McPt.GetSize();
          //std::cout<<"nDataParticles: "<<nDataParticles<<"   nMcParticles: "<<nMcParticles<<std::endl;
	  ///Avoid different final state comparison
	  if( nDataParticles != nMcParticles ) continue;

	  Double_t event_distance_Med= -99;
	  Double_t SumMed_dPt2 = 0, SumMed_dEta2 = 0, SumMed_dPhi2 = 0;

          //std::cout<<"Starting particles matching"<<std::endl;
	  for(int idt = 0; idt < nDataParticles; ++idt){
	    bool repaired = false;    
	    Int_t nsame_flavor = 0;
	    Double_t tmp_dPt = 0, tmp_dEta = 0, tmp_dPhi = 0;

	    for(Int_t imc = 0; imc < nMcParticles; ++imc){	      
	      ///Avoid different Data-MC particles comparison
	      if(SetFlavorConstraint == "true" && DataId[idt] != McId[imc]) continue;

	      ///Compute preliminary particles distance
	      Double_t dPt  = fabs(DataPt[idt]-McPt[imc])*ScaledPt;
	      Double_t dEta = fabs(DataEta[idt]-McEta[imc])*ScaledEta;
	      Double_t dPhi = fabs(DataPhi[idt]-McPhi[imc])*ScaledPhi;
              //std::cout<<"Deltas (Pt/Eta/Phi): "<<dPt<<"/"<<dEta<<"/"<<dPhi<<std::endl;
	      if(SetFlavorConstraint == "true"){
		if(nsame_flavor == 0){
		  tmp_dPt  = dPt;
		  tmp_dEta = dEta;
		  tmp_dPhi = dPhi;
		}
		else{
		  ///Repair the previous one
		  if(repaired != true){
		    tmp_dPt  = 0.5*tmp_dPt;
		    tmp_dEta = 0.5*tmp_dEta;
		    tmp_dPhi = 0.5*tmp_dPhi;
		    repaired = true;
		  }
		  ///Append the new one
		  tmp_dPt  += 0.5*dPt;
		  tmp_dEta += 0.5*dEta;
		  tmp_dPhi += 0.5*dPhi;
		}
		nsame_flavor++;
	      }
	      else{
		tmp_dPt  += 0.5*dPt;
		tmp_dEta += 0.5*dEta;
		tmp_dPhi += 0.5*dPhi;
	      }
	      
	      if( Verbose == 3 )
		std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<std::endl;
	    }///Ends Data event loop
	  
	    SumMed_dPt2  += tmp_dPt*tmp_dPt;
	    SumMed_dEta2 += tmp_dEta*tmp_dEta;
	    SumMed_dPhi2 += tmp_dPhi*tmp_dPhi;
	  }///Ends MC event loop
	
	  ///Compute final Data-MC events distance & searches for minimum distance
          //std::cout<<"SumMedPt2: "<<SumMed_dPt2<<"  SumMedEta2: "<<SumMed_dEta2<<"   SumMedPhi2: "<<SumMed_dPhi2<<std::endl;
	  if(SumMed_dPt2 > 0){
	    event_distance_Med = sqrt(SumMed_dPt2 + SumMed_dEta2 + SumMed_dPhi2);
	    if(event_distance_Med < min_distance_Med){
	      min_distance_Med = event_distance_Med;
	      imc_min = mc;
	      f_type = *McType;
	    }
	    if( Verbose > 2 ) std::cout<<"Event_distance (Mean) = "<<event_distance_Med<<std::endl;  
	  }
        }///End MC sample loop

        MinDistance.push_back(min_distance_Med);
        McFile.push_back(f_type);
        PairedMc.push_back(imc_min);
	if( Verbose > 1 ) std::cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Med<<std::endl;
      }///End Data sample loop
  
      fme_tree->Fill();
      delete fData;
    }///Ends the loop over data files
    

    t2.Stop();
    std::cout<<"Finished MSC..."<<std::endl;
    return fme_tree;
  };
  //===================================================================================================================================================================
  
  
  
  
  
  ///Spliting samples in mini-batches to comport them to the number of cores
  Int_t nMcSamples = vMCs.size();
  Int_t N_Cores = NCores;
  const Int_t PrimaryDivision = nMcSamples/N_Cores;
  Int_t Resting = nMcSamples % N_Cores;
  Int_t nBatches = (Resting >= 1)? PrimaryDivision+1 : PrimaryDivision;
  TList *list = new TList;  
  
  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Mapping Event to Event"<<ansi_blue<<" ]::::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
  for(Int_t ib = 0; ib < nBatches; ib++){
    std::cout<<"\n:: Processing MC Batch: "<<ib<<std::endl;
    std::vector<std::string> McBatches;

    if(ib < PrimaryDivision){
      for(Int_t iS = 0; iS < (Int_t)N_Cores; iS++)
        McBatches.push_back( vMCs[ib*N_Cores+iS] );
      TProcPool workers(N_Cores);
      if(PhSDrMethod == "mindr") list->Add( (TTree*)workers.ProcTree(McBatches, MPD) );
      if(PhSDrMethod == "mean")  list->Add( (TTree*)workers.ProcTree(McBatches, MSC) );
    }
    else{
      for(Int_t iS = 0; iS < (Int_t)Resting; iS++)
	McBatches.push_back( vMCs[ib*N_Cores+iS] );
      TProcPool workers(Resting);
      if(PhSDrMethod == "mindr") list->Add( (TTree*)workers.ProcTree(McBatches, MPD) );
      if(PhSDrMethod == "mean")  list->Add( (TTree*)workers.ProcTree(McBatches, MSC) );
    }
  }


  ///Merge the trees
  TTree *phs_tree = TTree::MergeTrees(list);
  phs_tree->SetName("CartographerResults");


  ///Cleaning the temporary folder
  gSystem->Exec("rm -r FME_USAGE");

  
  ///Send final tree merged from trees coming from all used cores
  return phs_tree;
}



#endif
