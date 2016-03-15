///####################################################################################################################
///=========================================[ Fast Matrix Element Module ]=============================================
///=====================================[ Code Author: Miqueias M. de Almeida ]========================================
///####################################################################################################################


#ifndef ComputePhsDR_h
#define ComputePhsDR_h

#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>

#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TStopwatch.h>

///Headers to TProcPool
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TProcPool.h>
#include <PoolUtils.h>
#include <TSystem.h>



void FindScaleFactors(FmeSetup Setup, Double_t *f_scale_dPt, Double_t *f_scale_dEta){
  TH1D *stackpt  = new TH1D("stackpt","",10000,0,10000);
  TH1D *stacketa = new TH1D("stacketa","",10000,0,10000);
  
  Double_t pt_sum = 0, eta_sum = 0, total = Setup.vMCs.size();
  for(Int_t isample=0; isample<(Int_t)Setup.vMCs.size(); isample++){
    TFile *ftmp = TFile::Open((TString)Setup.vMCs[isample]);
    TTree *ttmp = (TTree*)ftmp->Get(Setup.TTreeName);
    
    ttmp->Project("stackpt", Setup.PtBranch);
    if(Setup.ScaleMethod == "mean")    pt_sum += stackpt->GetMean();
    if(Setup.ScaleMethod == "maximum") pt_sum += stackpt->GetMaximum();
    ttmp->Project("stacketa", Setup.EtaBranch);
    if(Setup.ScaleMethod == "mean")    eta_sum += stacketa->GetMean();
    if(Setup.ScaleMethod == "maximum") eta_sum += stacketa->GetMaximum();
  }
  
  *f_scale_dPt  = pt_sum/total;
  *f_scale_dEta = eta_sum/total;
  
  std::cout<<Form(":: [NOTE] Setting scale_dPt = %.3f and scale_dEta = %.3f",*f_scale_dPt,*f_scale_dEta)<<std::endl;
  return;
}

///######################################## FastME Main Function ######################################################
TTree *ComputePhsDR(FmeSetup Setup){

  std::cout<<"::::::::::::::::::::::::::::::::[ Getting User Configuration ]:::::::::::::::::::::::::::::::::"<<std::endl;
  TFile *fData 			= Setup.DataFile;
  Int_t nData 			= Setup.NData;
  TString TreeName 		= Setup.TTreeName;
  TString McType_branch 	= Setup.McTypeBranch;
  TString Id_branch 		= Setup.IdBranch;
  TString Pt_branch 		= Setup.PtBranch;
  TString Eta_branch 		= Setup.EtaBranch;
  std::vector<std::string> MCs 	= Setup.vMCs;
  Int_t N_MCT			= Setup.NMCT;
  UInt_t N_Cores		= Setup.NCores;
  Int_t N_FSParticles		= Setup.NFSParticles;
  TString PhSDr_Method		= Setup.PhSDrMethod;
  TString FlavorConstraint	= Setup.SetFlavorConstraint;
  Float_t MC_Limit		= Setup.MCLimit;
  Double_t scale_dPt		= Setup.Scale_dPt;
  Double_t scale_dEta		= Setup.Scale_dEta;
  Int_t verbose			= Setup.Verbose;

  std::cout<<"Initials scale_dPt and scale_dEta -----> "<<scale_dPt<<", "<<scale_dEta<<std::endl;  
  if(scale_dPt == -1 || scale_dEta == -1) FindScaleFactors(Setup, &scale_dPt, &scale_dEta);

  ///Timming full process
  std::cout<<"::::::::::::::::::::::::::::::::[ Computing Events Distance ]:::::::::::::::::::::::::::::::::"<<std::endl;
  
  TStopwatch t1;
  t1.Start();
  
  ///TProcPool declaration to objects to be analised  
  auto workItem = [fData, nData, TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch, N_MCT, N_FSParticles,
		   PhSDr_Method, FlavorConstraint, MC_Limit, scale_dPt, scale_dEta, verbose]
		   (TTreeReader &tread) -> TObject* {
    TStopwatch t2;
        
    ///Addresses the MC branches to be used
    TTreeReaderValue<Int_t>    McType(tread, McType_branch); ///McType for Signal=0 and Background >0
    TTreeReaderArray<Int_t>    McId(tread, Id_branch);
    TTreeReaderArray<Double_t> McPt(tread, Pt_branch);
    TTreeReaderArray<Double_t> McEta(tread, Eta_branch);
    
    ///Addresses the Data branches to be used
    TTreeReader refReader(TreeName,fData);
    TTreeReaderArray<Int_t>    DataId(refReader, Id_branch);
    TTreeReaderArray<Double_t> DataPt(refReader, Pt_branch);
    TTreeReaderArray<Double_t> DataEta(refReader, Eta_branch);


    ///Tree to store the results from analysis
    Int_t iEvent, TMcType, Indice;
    Double_t Mdist;
    TTree *fme_tree = new TTree("fme_tree","temporary");
    fme_tree->SetDirectory(0);
    fme_tree->Branch("iEvent",&iEvent,"iEvent/I");
    fme_tree->Branch("Mdist",&Mdist,"Mdist/D");
    fme_tree->Branch("TMcType",&TMcType,"TMcType/I");
    fme_tree->Branch("Indice",&Indice,"Indice/I");
    
    ///Loop on Data events
    for(Int_t dt=0; dt<nData; dt++){
      if( verbose != 0 && ((dt!= 0 && dt%(nData/10) == 0) || (nData-dt) == 1) ){ 
	std::cout<< Form(":: [Remaining]:  %i Events\t\t[Elapsed]:  ",nData-dt);
	t2.Stop();
	t2.Print();
	t2.Continue();
      }
      refReader.SetEntry(dt); ///Move on Data loop
      Double_t min_distance_Min = 1.e15;
      Double_t min_distance_Med = 1.e15;
      Int_t imc_min = -1;
      Int_t f_type=-99;
      Int_t nMonteCarlo = tread.GetEntries(true);
      if(MC_Limit != -1 && MC_Limit >= 1) nMonteCarlo = MC_Limit;
      if(MC_Limit != -1 && MC_Limit < 1) nMonteCarlo = (Int_t)(MC_Limit*nMonteCarlo);
      
      for(Int_t mc=0; mc<nMonteCarlo; mc++){
	tread.SetEntry(mc); ///Move on MC loop
        bool acept = true;
	
      ///==============================================================================================================
      ///::::::::::::::::::::::::: Fast Matrix Element methods to compute Data - MC distance ::::::::::::::::::::::::::
      ///::::::::::::::::::::::::::::::::::::: Currently 2 Methods Available ::::::::::::::::::::::::::::::::::::::::::
      ///==============================================================================================================
	Double_t event_distance_Min= -99, event_distance_Med= -99;
	Double_t SumMed_dPt2 = 0, SumMed_dEta2 = 0;
	Double_t SumMin_dPt2 = 0, SumMin_dEta2 = 0;
	
	//stores flags to sinalize when a MC object is already selected
	std::vector<int> vmin_imc;
	for(int sl=0; sl<N_FSParticles; sl++) vmin_imc.push_back(-1);
	
	for(int idt=0; idt<N_FSParticles; idt++){
	  Double_t min_particles_distance = 1.E15;
	  Double_t particles_distance = -1.;
	  int min_imc = -1;
    
	  Int_t nsame_flavor = 0;
	  Double_t tmp_dPt = 0, tmp_dEta = 0;
	  for(int imc=0; imc<N_FSParticles; imc++){
	    ///Avoid different Data-MC particles comparison
	    if(FlavorConstraint == "true" && DataId[idt] != McId[imc]) continue;
	    ///Avoid leptons-jets comparison
	    else if(FlavorConstraint == "false"){
	      if(abs(DataId[idt])== 11 && (abs(McId[imc])!= 11 || abs(McId[imc])!= 13)) continue;
	      if(abs(DataId[idt])== 13 && (abs(McId[imc])!= 11 || abs(McId[imc])!= 13)) continue;
	      if(abs(McId[idt])== 11 && (abs(DataId[imc])!= 11 || abs(DataId[imc])!= 13)) continue;
	      if(abs(McId[idt])== 13 && (abs(DataId[imc])!= 11 || abs(DataId[imc])!= 13)) continue;
	    }
	    ///Compute preliminary particles distance
	    Double_t dPt  = (DataPt[idt]-McPt[imc])/(scale_dPt);
	    Double_t dEta = (DataEta[idt]-McEta[imc])/(scale_dEta);
	

	  ///_______________________ For proximity comparison method __________________________________________________
	    if( PhSDr_Method == "mindr")
	      if(vmin_imc[imc] == -1){
		particles_distance = sqrt(dPt*dPt + dEta*dEta);
		if( verbose == 3 ) std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<"   part_dist: "<<particles_distance<<std::endl;
		if(particles_distance < min_particles_distance){
		  min_imc = imc;
		  min_particles_distance = particles_distance;
		}
	      }
	  ///¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
	  
	  
  	  ///________________________ Only for media comparison method ________________________________________________
	    if(PhSDr_Method == "media"){
	      if(nsame_flavor == 0){
		tmp_dPt  = dPt/scale_dPt;
		tmp_dEta = dEta/scale_dEta;
	      }
	      else{
		///Repair the previous one
		SumMed_dPt2  += pow(0.5*tmp_dPt,2);
		SumMed_dEta2 += pow(0.5*tmp_dEta,2);
		///Append the new one
		SumMed_dPt2  += pow(0.5*tmp_dPt,2);
		SumMed_dEta2 += pow(0.5*tmp_dEta,2);
	      }
	      nsame_flavor++;
	      if( verbose == 3 )
		std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<std::endl;
	    }
	  ///¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨

	  }///Ends MC event loop
	  
	  if(PhSDr_Method == "mindr"){
	    ///Monitor of chosen MCs to avoid object recounting and wrong pairing
	    if(min_imc == -1){
	      acept = false;
	      break;
	    }
	    vmin_imc[min_imc] = 1;//changes the flag for current MC object
	    if( verbose == 3 ) std::cout<<"Chosen->>  MCPos: "<<min_imc<<"   ID: "<<McId[min_imc]<<std::endl;
	    ///For proximity comparison method
	    SumMin_dPt2 += pow( (DataPt[idt]-McPt[min_imc])/(scale_dPt), 2 );
	    SumMin_dEta2 += pow( (DataEta[idt]-McEta[min_imc])/(scale_dEta), 2 );
	  }
	  
	  if(PhSDr_Method == "media" && nsame_flavor == 1){
	    SumMed_dPt2  += tmp_dPt*tmp_dPt;
	    SumMed_dEta2 += tmp_dEta*tmp_dEta;
	  }
	}///Ends Data event loop
	
	///Compute final Data-MC events distance & searches for minimum distance
	if(PhSDr_Method == "mindr" && acept == true){
	  event_distance_Min = sqrt(SumMin_dPt2 + SumMin_dEta2);
	  if(event_distance_Min < min_distance_Min){
	    min_distance_Min = event_distance_Min;
	    imc_min = mc;
	  }
	  if( verbose > 2 ) std::cout<<"Event_distance(MinDr) = "<<event_distance_Min<<std::endl;  
	}
	
	if(PhSDr_Method == "media" && SumMed_dPt2 > 0){
	  event_distance_Med = sqrt(SumMed_dPt2 + SumMed_dEta2);
	  if(event_distance_Med < min_distance_Med){
	    min_distance_Med = event_distance_Med;
	    imc_min = mc;
	  }
	  if( verbose > 2 ) std::cout<<"Event_distance (Media) = "<<event_distance_Med<<std::endl;  
	}
	///Stores the MC type
	f_type = *McType;
      ///================================================================================================================
      ///================================================================================================================
	
      }///End MC sample loop

      
      ///Stores the minimum distances found
      iEvent = dt;
      if(PhSDr_Method == "mindr"){
	Mdist = min_distance_Min;
	TMcType = f_type;
	Indice  = imc_min;
	if( verbose > 1 ) std::cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Min<<std::endl;
      }
      if(PhSDr_Method == "media"){
	Mdist = min_distance_Med;
	TMcType = f_type;
	Indice  = imc_min;
	if( verbose > 1 ) std::cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Med<<std::endl;
      }
      
      fme_tree->Fill();
    }///End Data sample loop
    
    t2.Stop();
    delete fData;
    return fme_tree;
  };
  
  ///Calls analysis through TProcPool
  TProcPool workers(N_Cores);
  TTree *mtree = (TTree*)workers.ProcTree(MCs, workItem);

  ///________________________________ Stoping timming ________________________________________________________
  std::cout<<"\n::::::::::::::::::::::::::::::::::::[ Process Finished ]::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<":: [Analysis Total Time]: "; t1.Stop(); t1.Print();
  std::cout<<":: [Sending TTree results...]"<<std::endl;
  std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"<<std::endl;
  ///---------------------------------------------------------------------------------------------------------

  ///Final tree merged from trees coming from all cores used
  return mtree;
}

#endif
