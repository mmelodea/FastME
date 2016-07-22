///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::[ PhsDrComputer - Compute Events Distance ]::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef ComputePhsDR_h
#define ComputePhsDR_h

#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <cmath>

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



void FindScaleFactors(FmeSetup Setup, Double_t *f_scale_dPt, Double_t *f_scale_dEta){
  gROOT->SetBatch();
  TCanvas *temp = new TCanvas();
  Double_t pt_sum = 0, eta_sum = 0, total = Setup.vMCs.size();
  for(Int_t isample=0; isample<(Int_t)total; isample++){
    TFile *ftmp = TFile::Open((TString)Setup.vMCs.at(isample));
    TTree *ttmp = (TTree*)ftmp->Get(Setup.TTreeName);
    
    if(*f_scale_dPt < 0){
      TString draw_pt = Setup.PtBranch+" >> stackpt";
      ttmp->Draw(draw_pt);
      TH1D *stackpt = (TH1D*)gDirectory->Get("stackpt");
      if(*f_scale_dPt == -1) pt_sum += stackpt->GetMean();
      if(*f_scale_dPt == -2) pt_sum += stackpt->GetBinCenter(stackpt->GetMaximumBin());
    }
    if(*f_scale_dEta < 0){
      TString draw_eta = Setup.EtaBranch+" >> stacketa";
      ttmp->Draw(draw_eta);
      TH1D *stacketa = (TH1D*)gDirectory->Get("stacketa");
      if(*f_scale_dEta == -1) eta_sum += fabs(stacketa->GetMean());//Only in case you are in region shifted from 0!
      if(*f_scale_dEta == -2) eta_sum += fabs(stacketa->GetBinCenter(stacketa->GetMinimum()));
    }
  }
  
  if(*f_scale_dPt < 0)  *f_scale_dPt  = pt_sum/total;
  if(*f_scale_dEta < 0) *f_scale_dEta = eta_sum/total;
  
  std::cout<<":: ["<<ansi_yellow<<"NOTE"<<ansi_reset<<Form("] Setting scale_dPt = %.3f and scale_dEta = %.3f",*f_scale_dPt,*f_scale_dEta)<<std::endl;
  temp->Close();
  return;
}





///######################################## FastME Main Function ######################################################
TTree *ComputePhsDR(FmeSetup Setup){

  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Getting User Configuration"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
  std::vector<std::string> Datas = Setup.vDatas;
  Int_t nData 			 = Setup.DTLimit;
  TString TreeName 		 = Setup.TTreeName;
  TString McType_branch 	 = Setup.McTypeBranch;
  TString Id_branch 		 = Setup.IdBranch;
  TString Pt_branch 		 = Setup.PtBranch;
  TString Eta_branch 		 = Setup.EtaBranch;
  std::vector<std::string> MCs 	 = Setup.vMCs;
  Int_t N_MCT			 = Setup.NMCT;
  UInt_t N_Cores		 = Setup.NCores;
  TString PhSDr_Method		 = Setup.PhSDrMethod;
  TString FlavorConstraint	 = Setup.SetFlavorConstraint;
  Float_t MC_Limit		 = Setup.MCLimit;
  Double_t scale_dPt		 = Setup.ScaledPt;
  Double_t scale_dEta		 = Setup.ScaledEta;
  Int_t verbose			 = Setup.Verbose;


  ///Verifying the scale factors
  if(scale_dPt < 0 || scale_dEta < 0){
    std::cout<<":: ["<<ansi_yellow<<"Initials scale_dPt and scale_dEta -----> "<<scale_dPt<<", "<<scale_dEta<<" @@Computing new scale factors..."<<ansi_reset<<"]"<<std::endl;
    FindScaleFactors(Setup, &scale_dPt, &scale_dEta);
  }

  
  ///Timming full process
  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Computing Events Distance"<<ansi_blue<<" ]::::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;
  
  
  ///TProcPool declaration to objects to be analised  
  auto workItem = [Datas, nData, TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch, 
		   N_MCT, PhSDr_Method, FlavorConstraint, MC_Limit, scale_dPt, scale_dEta, verbose]
		   (TTreeReader &tread) -> TObject* {
		     
    std::cout<<ansi_yellow<<"::----->>> Activating core <<<-----::"<<ansi_reset<<std::endl;
    TStopwatch t2;
        
    ///Addresses the MC branches to be used
    TTreeReaderValue<int>    	McType(tread, McType_branch); ///McType for Signal=0 and Background >0
    TTreeReaderArray<int>	McId(tread, Id_branch);
    TTreeReaderArray<double>	McPt(tread, Pt_branch);
    TTreeReaderArray<double>	McEta(tread, Eta_branch);


    ///Tree to store the results from analysis
    Int_t DataFile;
    std::vector<double> Mdist;
    std::vector<int> Indice, TMcType, DtObjFlag, fDtObjFlag; 
    TTree *fme_tree = new TTree("fme_tree","From FastME Phase Space Analysis");
    fme_tree->SetDirectory(0);
    fme_tree->Branch("DataFile",&DataFile,"DataFile/I");
    fme_tree->Branch("MinDistance","std::vector<double>",&Mdist);
    fme_tree->Branch("PairedMC","std::vector<int>",&Indice);
    fme_tree->Branch("PairedMCType","std::vector<int>",&TMcType);
    fme_tree->Branch("DataObjFlag","std::vector<int>",&fDtObjFlag);


    ///Loop over the different data files
    for(Int_t idata=0; idata<(Int_t)Datas.size(); idata++){
      DataFile = idata;
      Mdist.clear();
      Indice.clear();
      TMcType.clear();

      TFile *fData = TFile::Open( (TString)Datas.at(idata) );
      TTreeReader refReader(TreeName,fData);
      TTreeReaderArray<int>      DataId(refReader, Id_branch);
      TTreeReaderArray<double>   DataPt(refReader, Pt_branch);
      TTreeReaderArray<double>   DataEta(refReader, Eta_branch);

    
      ///Loop on Data events
      Int_t nDataEv = refReader.GetEntries(true);
      if(nData != -1 && nData >= 1 && nData < nDataEv) nDataEv = nData;
      for(Int_t dt=0; dt<nDataEv; dt++){
	refReader.SetEntry(dt); ///Move on Data loop                                                              
      
	if( verbose != 0 && ((dt!= 0 && nDataEv > 10 && dt%(nDataEv/10) == 0) || (nDataEv-dt) == 1) ){
	  std::cout<<":: ["<<ansi_violet<<"Remaining from "<<fData->GetName()<<" for MC "<<*McType<<"/Elapsed"<<ansi_reset<<"]:    "<<nDataEv-dt<<"/ "<<(Int_t)t2.RealTime()<<"seg";
	  t2.Continue();
          if(nDataEv-dt == 1) std::cout<<ansi_yellow<<" <<-(Concluded!)"<<ansi_reset<<std::endl;
	  else std::cout<<std::endl;
	}

	Double_t min_distance_Min = 1.e15;
	Double_t min_distance_Med = 1.e15;
	Int_t imc_min = -1;
	Int_t f_type=-99;
	Int_t nMonteCarlo = tread.GetEntries(true);
	if(MC_Limit != -1 && MC_Limit >= 1) nMonteCarlo = MC_Limit;
	if(MC_Limit != -1 && MC_Limit < 1) nMonteCarlo = (Int_t)(MC_Limit*nMonteCarlo);

	for(Int_t mc=0; mc<nMonteCarlo; mc++){
	  //Initialyze the vector
	  DtObjFlag.clear();
	  for(int sl=0; sl<(int)DataId.GetSize(); sl++) DtObjFlag.push_back(-1);

	  tread.SetEntry(mc); ///Move on MC loop
	  bool acept = true;
	
	  ///==============================================================================================================
	  ///::::::::::::::::::::::::: Fast Matrix Element methods to compute Data - MC distance ::::::::::::::::::::::::::
	  ///::::::::::::::::::::::::::::::::::::: Currently 2 Methods Available ::::::::::::::::::::::::::::::::::::::::::
	  ///==============================================================================================================
	  Double_t event_distance_Min= -99, event_distance_Med= -99;
	  Double_t SumMed_dPt2 = 0, SumMed_dEta2 = 0;
	  Double_t SumMin_dPt2 = 0, SumMin_dEta2 = 0;
	
	  for(int imc=0; imc<(int)McId.GetSize(); imc++){
	    Double_t min_particles_distance = 1.E15;
	    Double_t particles_distance = -1.;
	    int sel_data = -1;

	    bool repaired = false;    
	    Int_t nsame_flavor = 0;
	    Double_t tmp_dPt = 0, tmp_dEta = 0;
	    for(int idt=0; idt<(int)DataId.GetSize(); idt++){
	      if(PhSDr_Method == "mindr" && DtObjFlag.at(idt) == 1) continue;///Skip data object already selected
	      
	      ///Avoid different Data-MC particles comparison
	      if(FlavorConstraint == "true" && DataId[idt] != McId[imc]) continue;
	      ///Avoid leptons-jets comparison
	      else if(FlavorConstraint == "false"){
		if( (abs(DataId[idt])!= 11 && abs(DataId[idt])!= 13) && (abs(McId[imc])== 11 || abs(McId[imc])== 13) ) continue;
		if( (abs(DataId[idt])== 11 || abs(DataId[idt])== 13) && (abs(McId[imc])!= 11 && abs(McId[imc])!= 13) ) continue;
	      }


	      ///Compute preliminary particles distance
	      Double_t dPt  = (DataPt[idt]-McPt[imc])/(scale_dPt);
	      Double_t dEta = (DataEta[idt]-McEta[imc])/(scale_dEta);
	

	      ///_______________________ For proximity comparison method __________________________________________________
	      if( PhSDr_Method == "mindr"){
		particles_distance = sqrt(dPt*dPt + dEta*dEta);
		if( verbose == 3 ) std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<"   part_dist: "<<particles_distance<<std::endl;
		if(particles_distance < min_particles_distance){
		  sel_data = idt;
		  min_particles_distance = particles_distance;
		}
	      }
	      ///----------------------------------------------------------------------------------------------------------
	  
	  
	      ///________________________ Only for mean comparison method ________________________________________________
	      if(PhSDr_Method == "mean"){
		if(FlavorConstraint == "true"){
		  if(nsame_flavor == 0){
		    tmp_dPt  = dPt;
		    tmp_dEta = dEta;
		  }
		  else{
		    ///Repair the previous one
		    if(repaired != true){
		      tmp_dPt  = 0.5*tmp_dPt;
		      tmp_dEta = 0.5*tmp_dEta;
		      repaired = true;
		    }
		    ///Append the new one
		    tmp_dPt  += 0.5*dPt;
		    tmp_dEta += 0.5*dEta;
		  }
		  nsame_flavor++;
		}
		else{
		  tmp_dPt  += 0.5*dPt;
		  tmp_dEta += 0.5*dEta;
		}
		if( verbose == 3 )
		  std::cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<std::endl;
	      }
	      ///--------------------------------------------------------------------------------------------------------

	    }///Ends Data event loop
	  
	    if(PhSDr_Method == "mindr"){
	      ///Monitor of chosen MCs to avoid object recounting and wrong pairing
	      if(sel_data == -1){
		acept = false;
		break;
	      }
	      DtObjFlag[sel_data] = 1;///changes the flag for current Data object

	      if( verbose == 3 ) std::cout<<"Chosen->>  DtPos: "<<sel_data<<"   ID: "<<DataId[sel_data]<<std::endl;
	      ///For proximity comparison method
	      SumMin_dPt2 += pow( (DataPt[sel_data]-McPt[imc])/(scale_dPt), 2 );
	      SumMin_dEta2 += pow( (DataEta[sel_data]-McEta[imc])/(scale_dEta), 2 );
	    }
	  
	    if(PhSDr_Method == "mean"){
	      SumMed_dPt2  += tmp_dPt*tmp_dPt;
	      SumMed_dEta2 += tmp_dEta*tmp_dEta;
	    }
	  }///Ends MC event loop
	
	  ///Compute final Data-MC events distance & searches for minimum distance
	  if(PhSDr_Method == "mindr" && acept == true){
	    event_distance_Min = sqrt(SumMin_dPt2 + SumMin_dEta2);
	    if(event_distance_Min < min_distance_Min){
	      min_distance_Min = event_distance_Min;
	      imc_min = mc;
	      fDtObjFlag = DtObjFlag;
	    }
	    if( verbose > 2 ) std::cout<<"Event_distance (MinDr) = "<<event_distance_Min<<std::endl;  
	  }
	
	  if(PhSDr_Method == "mean" && SumMed_dPt2 > 0){
	    event_distance_Med = sqrt(SumMed_dPt2 + SumMed_dEta2);
	    if(event_distance_Med < min_distance_Med){
	      min_distance_Med = event_distance_Med;
	      imc_min = mc;
	    }
	    if( verbose > 2 ) std::cout<<"Event_distance (Mean) = "<<event_distance_Med<<std::endl;  
	  }

	  ///Stores the MC type
	  f_type = *McType;
	  ///================================================================================================================
	  ///================================================================================================================
	
      }///End MC sample loop

      
	///Stores the minimum distances found
	if(PhSDr_Method == "mindr"){
	  Mdist.push_back(min_distance_Min);
	  TMcType.push_back(f_type);
	  Indice.push_back(imc_min);
	  if( verbose > 1 ) std::cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Min<<std::endl;
	}
	if(PhSDr_Method == "mean"){
          Mdist.push_back(min_distance_Med);
          TMcType.push_back(f_type);
          Indice.push_back(imc_min);
	  if( verbose > 1 ) std::cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Med<<std::endl;
	}
	
	//fme_tree->Fill();
      }///End Data sample loop
  
      fme_tree->Fill();
      delete fData;
    }///Ends the loop over data files
    

    t2.Stop();
    return fme_tree;
  };

  


  ///Calls analysis through TProcPool
  TProcPool workers(N_Cores);
  TTree *mtree = (TTree*)workers.ProcTree(MCs, workItem);
  

  ///________________________________ Stoping timming ________________________________________________________
  std::cout<<ansi_blue<<std::endl;
  std::cout<<":::::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Process Finished"<<ansi_blue<<" ]::::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<":: ["<<ansi_cyan<<"Sending PhsDrComputer Results to Discriminant Computer"<<ansi_blue<<"]"<<std::endl;
  std::cout<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<ansi_reset<<std::endl;
  ///---------------------------------------------------------------------------------------------------------

  
  ///Send final tree merged from trees coming from all cores used
  return mtree;
}

#endif
