///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::[ Librarian - Handles the MC files indexing ]:::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef Librarian_h
#define Librarian_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>


void Indexer(FmeSetup *Setup){

  ///Create a folder to keep data to be used by FastME
  ///Only the important branches and the defined events quantity are keept
  gSystem->Exec("mkdir FME_USAGE");
  gSystem->Exec("mkdir FME_USAGE/DATA");
  gSystem->Exec("mkdir FME_USAGE/MC_TEMPLATES");


  ///Preparing the MC templates
  std::cout<<":: Reducing data files: "<<std::endl;  
  std::cout<<"::------------------------------------------------------------::"<<std::endl;
  Int_t nDataFiles = Setup->vDatas.size();
  for(Int_t ifile = 0; ifile < nDataFiles; ++ifile){

    TFile *org_file = TFile::Open( (TString)Setup->vDatas[ifile] );
    TTree *org_tree = (TTree*)org_file->Get(Setup->TTreeName);
    org_tree->SetBranchStatus("*",kFALSE);
    org_tree->SetBranchStatus(Setup->IdBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->PtBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->EtaBranch,kTRUE);

    TFile *fin_file = new TFile( Form("FME_USAGE/DATA/Reduced_file_from_original_data_file_%i.root",ifile), "recreate" );
    TTree *fin_tree = org_tree->CloneTree();
    org_file->Close();
    
    fin_tree->Write();
    fin_file->Close();
    std::cout<<ansi_green<<"[REDUCED] "<<ansi_reset<<Setup->vDatas[ifile]<<std::endl;
    Setup->vDatas[ifile] = Form("FME_USAGE/DATA/Reduced_file_from_original_data_file_%i.root",ifile);
  }

    


  
  
  ///Preparing the MC templates
  std::cout<<":: Indexing MC files: "<<std::endl;
  std::cout<<"::------------------------------------------------------------::"<<std::endl;
  Int_t nMcFiles = Setup->vMCs.size();
  for(Int_t ifile = 0; ifile < nMcFiles; ++ifile){

    TFile *org_file = TFile::Open( (TString)Setup->vMCs[ifile] );
    TTree *org_tree = (TTree*)org_file->Get(Setup->TTreeName);
    org_tree->SetBranchStatus("*",kFALSE);
    org_tree->SetBranchStatus(Setup->IdBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->PtBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->EtaBranch,kTRUE);

    TFile *fin_file = new TFile( Form("FME_USAGE/MC_TEMPLATES/Indexed_file_from_original_MC_file_%i.root",ifile), "recreate" );
    TTree *fin_tree = org_tree->CloneTree();
    org_file->Close();

    Int_t McFileIndex = ifile;
    TBranch *file_index = fin_tree->Branch("McFileIndex",&McFileIndex,"McFileIndex/I");
    Int_t Nentries = fin_tree->GetEntries();


    for(Int_t i = 0; i < Nentries; ++i) file_index->Fill();


    fin_tree->Write();
    fin_file->Close();
    std::cout<<ansi_green<<"[INDEXED] "<<ansi_reset<<Setup->vMCs[ifile]<<std::endl;
    Setup->vMCs[ifile] = Form("FME_USAGE/MC_TEMPLATES/Indexed_file_from_original_MC_file_%i.root",ifile);
  }




  return;

}



#endif
