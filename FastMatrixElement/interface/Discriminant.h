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
  

  Int_t nDtFiles	= Setup.vDatas.size(); //Number of data files
  Int_t nMcFiles	= Setup.vMCs.size(); //Number of MC files
  Int_t verbose		= Setup.Verbose;  

  TStopwatch t2;
  t2.Start();
  std::cout<<ansi_blue<<":::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Computing Discriminant"<<ansi_blue<<" ]::::::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;

  ///Set the input tree
  Int_t DtFile;
  std::vector<int> *TMcType=0, *Indice=0;
  std::vector<double> *Mdist=0;
  std::vector<int> *oDtObjFlag=0;
  mtree->SetBranchAddress("DataFile",&DtFile);
  mtree->SetBranchAddress("MinDistance",&Mdist);
  mtree->SetBranchAddress("PairedMCType",&TMcType);
  mtree->SetBranchAddress("PairedMC",&Indice);
  mtree->SetBranchAddress("DataObjFlag",&oDtObjFlag);  
  Int_t fentries = mtree->GetEntries();

  std::vector<Int_t> McIndex, McCat, DtObjFlag;
  Double_t Global_PsbDist;
  std::vector<Double_t> MinDist, Local_PsbDist;
  TTree *ftree = new TTree("FastME","Fast Matrix Element Analysis Results");
  ftree->SetDirectory(0);
  ftree->Branch("DataFile",&DtFile,"DataFile/I");
  ftree->Branch("PairedMC","std::vector<Int_t>",&McIndex);
  ftree->Branch("PairedMCType","std::vector<Int_t>",&McCat);
  ftree->Branch("DataObjFlag","std::vector<Int_t>",&DtObjFlag);
  ftree->Branch("MinDistance","std::vector<Double_t>",&MinDist);
  ftree->Branch("Global_PsbDist",&Global_PsbDist,"Global_PsbDist/D");
  ftree->Branch("Local_PsbDist","std::vector<Double_t> Local_PsbDist",&Local_PsbDist);

  ///Find the tree sectors (sections of tree for each MC)
  ///The data repeats for every MC type inserted in the analysis
  Int_t TreeSectors[nMcFiles];
  if(fentries % nMcFiles != 0){
    std::cout<<ansi_red<<"[Error]"<<ansi_reset<<" Something went wrong... non-integer tree sectors!!"<<std::endl;
    std::cout<<"Entries on internal tree: "<<fentries<<std::endl;
    std::cout<<"Number of MCs: "<<nMcFiles<<std::endl;
    std::cout<<"(fentries%nMcFiles): "<< fentries % nMcFiles <<std::endl;
    throw std::exception();
  }
                     
  for(Int_t ic=0; ic<(Int_t)nMcFiles; ic++) TreeSectors[ic] = ic*nDtFiles;
  
  ///Getting results from analysis
  ///Each tree row has all events from a data file
  std::cout<<":: Detected "<<nDtFiles<<" data files..."<<std::endl;
  for(Int_t ifile=0; ifile<nDtFiles; ifile++){


    ///Goes over the events from current file
    mtree->GetEntry(ifile);
    const Int_t nevents = Mdist->size();
    for(Int_t ievent=0; ievent<nevents; ievent++){
      if(ievent >= 10 && ievent%(nevents/10) == 0){
        std::cout<<":: ["<<ansi_violet<<"Remaning DataFile/Events"<<ansi_reset<<"]:  "<<Form("%i/%i/ %.3fseg",ifile,nevents-ievent,t2.RealTime())<<std::endl;
        t2.Continue();
      }
    
      Int_t MinSigIndex = -99, smin_ic = -1, bmin_ic = -1;
      Double_t min_dr_sig = 1.e15, global_min_dr_bkg = 1.e15;
      //Double_t local_min_dr_bkg[nMcFiles], local_min_bkg_index[nMcFiles];
      McIndex.clear();
      McCat.clear();
      MinDist.clear();
      Local_PsbDist.clear();
      DtObjFlag.clear();


      ///Looping over the MC sectors
      for(Int_t ic=0; ic<(Int_t)nMcFiles; ic++){
        mtree->GetEntry(ifile + TreeSectors[ic]); //TreeSectors aligns the results from different MCs

	if(verbose > 1)
          std::cout<<"Loading entry "<<ifile + TreeSectors[ic]<<"\tDataFile/Event/MCCat "<<ifile<<"/"<<ievent<<"/"<<(*TMcType).at(ievent)<<std::endl;
      
      
        ///Finds closest MC Signal
        if( (*TMcType).at(ievent) == 0 ){
   	  if( (*Mdist).at(ievent) < min_dr_sig ){
	    min_dr_sig = (*Mdist).at(ievent);
            MinDist.push_back( (*Mdist).at(ievent) );
	    MinSigIndex = (*Indice).at(ievent);
	    McCat.push_back( 0 );
	    smin_ic = ic;
	  }
        }

        ///Finds closest MC Background
        if( (*TMcType).at(ievent) > 0 ){
	  ///The general closest Background MC
          if( (*Mdist).at(ievent) < global_min_dr_bkg ){
            global_min_dr_bkg = (*Mdist).at(ievent);
            bmin_ic = ic;
	  }
          MinDist.push_back( (*Mdist).at(ievent) );
	  McCat.push_back( (*TMcType).at(ievent) );
	  ///Each MC Background
	  //if( (*Mdist).at(ievent) < local_min_dr_bkg[(*TMcType).at(ievent)] ){
	    //McCat[(*TMcType).at(ievent)] = (*TMcType).at(ievent);
	    //local_min_dr_bkg[(*TMcType).at(ievent)] = (*Mdist).at(ievent);
	    //local_min_bkg_index[(*TMcType).at(ievent)] = (*Indice).at(ievent);
	  //}
        }

      }//Ending the full verification for a event

      Global_PsbDist = GetPsbD(min_dr_sig, global_min_dr_bkg);

      //for(Int_t im=1; im<nMcFiles; im++){
        
        //Local_PsbDist[im] = GetPsbD(min_dr_sig, local_min_dr_bkg[im]);
      //}
      if( verbose > 1 )
        std::cout<< Form("GSigMin:   %f\t\tGBkgMin:   %f\t\tGPsbDMinDist:   %f", min_dr_sig, global_min_dr_bkg, Global_PsbDist) << std::endl;

      ///Set the final flags to the Data objects
      Int_t loadTTree_ic = (min_dr_sig < global_min_dr_bkg)? smin_ic : bmin_ic;
      mtree->GetEntry(ifile + TreeSectors[loadTTree_ic]);
      for(Int_t iobj=0; iobj<(Int_t)oDtObjFlag->size(); iobj++)
        DtObjFlag.push_back( (*oDtObjFlag).at(iobj) );

    
      ftree->Fill();
    }
  }
  
  ///________________________________ Stoping timming ________________________________________________________
  std::cout<<ansi_blue<<std::endl;
  std::cout<<"::::::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Process Finished"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<":: ["<<ansi_cyan<<"Computing Total Time"<<ansi_blue<<"]: "<<Form("%.3f seg", t2.RealTime())<<std::endl;
  std::cout<<":: ["<<ansi_cyan<<"Sending Discriminant Results"<<ansi_blue<<"]"<<std::endl;
  std::cout<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<ansi_reset<<std::endl;
  ///---------------------------------------------------------------------------------------------------------


  ///Clean memory from used pointers
  delete TMcType;
  delete Indice;
  delete Mdist;
  delete oDtObjFlag;



  ///Send the final tree to be stored with the full results from FastME analysis
  return ftree;

}


#endif
