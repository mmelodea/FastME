///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::[ Discriminant Module - Analyze PhsDrComputer's Output ]:::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Discriminant_h
#define Discriminant_h

#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <vector>

#include <TTree.h>
#include <TStopwatch.h>


///================ Discriminant Based on Distance ===================
Double_t GetPsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}
///===================================================================


///_______________________ Compute discriminant from MDMCED _________________________________________________
TTree *Discriminant(TTree *mtree, FmeSetup Setup){
  
  UInt_t N_Cores	= Setup.NCores;
  Int_t nData		= Setup.NData;
  const Int_t N_MCT	= Setup.NMCT;
  Int_t verbose		= Setup.Verbose;  

  TStopwatch t2;
  t2.Start();
  std::cout<<ansi_blue<<"::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Computing Discriminant"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;

  ///Set the input tree
  Int_t iEvent, TMcType, Indice;
  Double_t Mdist;
  mtree->SetBranchAddress("iEvent",&iEvent);
  mtree->SetBranchAddress("Mdist",&Mdist);
  mtree->SetBranchAddress("TMcType",&TMcType);
  mtree->SetBranchAddress("Indice",&Indice);  
  Int_t fentries = mtree->GetEntries();

  std::vector<Int_t> McIndex, McCat;
  Double_t Global_PsbDist;
  std::vector<Double_t> MinDist, Local_PsbDist;
  TTree *ftree = new TTree("FastME","Fast Matrix Element Analysis Results");
  ftree->SetDirectory(0);
  ftree->Branch("McIndex","std::vector<Int_t>",&McIndex);
  ftree->Branch("McCat","std::vector<Int_t>",&McCat);
  ftree->Branch("MinDist","std::vector<Double_t>",&MinDist);
  ftree->Branch("Global_PsbDist",&Global_PsbDist,"Global_PsbDist/D");
  ftree->Branch("Local_PsbDist","std::vector<Double_t> Local_PsbDist",&Local_PsbDist);

  ///Find the tree sectors
  Int_t TreeSectors[N_Cores];
  if(fentries % N_Cores != 0){
    std::cout<<ansi_red<<"[Error]"<<ansi_reset<<" Something gone wrong... non-integer tree sectors!!"<<std::endl;
    throw std::exception();
  }
  
  Int_t EndSector = fentries/N_Cores; //How many entries in each core job                                     
  for(Int_t ic=0; ic<(Int_t)N_Cores; ic++) TreeSectors[ic] = ic*EndSector;
  
  ///Getting results from analysis
  for(Int_t data=0; data<nData; data++){
    if( verbose != 0 && nData > 10 && data%(nData/10) == 0){
      std::cout<<":: ["<<ansi_violet<<"Remaining events"<<ansi_reset<<"]:  "<<nData-data<<"\t\t";
      t2.Print();
    }

    Int_t MinSigIndex = -99;// GMinBkgIndex = -99;
    Double_t min_dr_sig = 1.e15, global_min_dr_bkg = 1.e15;
    Double_t local_min_dr_bkg[N_MCT], local_min_bkg_index[N_MCT];
    McIndex.clear();
    McCat.clear();
    MinDist.clear();
    Local_PsbDist.clear();
    for(Int_t p=0; p<N_MCT; p++){
      local_min_bkg_index[p] = -99;
      local_min_dr_bkg[p] = 1.e15;
      McCat.push_back( -99 );
      McIndex.push_back( -99 );
      MinDist.push_back( -99. );
      Local_PsbDist.push_back( -99. );
    }
    for(Int_t ic=0; ic<(Int_t)N_Cores; ic++){
      mtree->GetEntry(TreeSectors[ic]+data);//TreeSectors aligns the results from different cores
      if(iEvent != data){
	std::cout<<ansi_red<<"[Error]"<<ansi_reset<<" Something gone wrong... iEvent != data"<<std::endl;
	throw std::exception();
      }
      
      ///Finds closet MC Signal
      if(TMcType == 0)
	if( Mdist < min_dr_sig ){
	  min_dr_sig = Mdist;
	  MinSigIndex = Indice;
	  McCat[0] = 0;
	}

      ///Finds closet MC Background
      if(TMcType > 0){
	///The general most close MC Background
        if( Mdist < global_min_dr_bkg ) global_min_dr_bkg = Mdist;
	///Each MC Background
	if( Mdist < local_min_dr_bkg[TMcType] ){
	  McCat[TMcType] = TMcType;
	  local_min_dr_bkg[TMcType] = Mdist;
	  local_min_bkg_index[TMcType] = Indice;
	}
      }
    }

    MinDist[0] = min_dr_sig;
    McIndex[0] = MinSigIndex;
    Global_PsbDist = GetPsbD(min_dr_sig, global_min_dr_bkg);

    for(Int_t im=1; im<N_MCT; im++){
      McIndex[im] = local_min_bkg_index[im];
      MinDist[im] = local_min_dr_bkg[im];
      Local_PsbDist[im] = GetPsbD(min_dr_sig, local_min_dr_bkg[im]);
    }
    if( verbose > 1 )
      std::cout<< Form("GSigMin:   %f\t\tGBkgMin:   %f\t\tGPsbDMinDist:   %f", min_dr_sig, global_min_dr_bkg, Global_PsbDist) << std::endl;
    
    ftree->Fill();
  }
  
  ///________________________________ Stoping timming ________________________________________________________
  std::cout<<ansi_blue<<std::endl;
  std::cout<<":::::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Process Finished"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<":: ["<<ansi_cyan<<"Computing Total Time"<<ansi_blue<<"]: "; t2.Stop(); t2.Print();
  std::cout<<":: ["<<ansi_cyan<<"Sending Discriminant Results"<<ansi_blue<<"]"<<std::endl;
  std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<ansi_reset<<std::endl;
  ///---------------------------------------------------------------------------------------------------------


  ///Send the final tree to be stored with the full results from FastME analysis
  return ftree;
}


#endif