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

  
  std::cout<<":: Indexing files: "<<std::endl;  
  gSystem->Exec("mkdir -p IndexedFiles");
  gSystem->Exec("rm -r IndexedFiles/");
  gSystem->Exec("mkdir IndexedFiles");

  for(Int_t ifile = 0; ifile < (Int_t)Setup->vMCs.size(); ifile++){

    TFile *org_file = TFile::Open( (TString)Setup->vMCs.at(ifile) );
    TTree *org_tree = (TTree*)org_file->Get(Setup->TTreeName);

    TFile *fin_file = new TFile( Form("IndexedFiles/Indexed_file_from_original_file_%i.root",ifile), "recreate" );
    TTree *fin_tree = org_tree->CloneTree();
    org_file->Close();

    Int_t McFileIndex = ifile;
    TBranch *file_index = fin_tree->Branch("McFileIndex",&McFileIndex,"McFileIndex/I");
    Int_t Nentries = fin_tree->GetEntries();


    for(Int_t i = 0; i < Nentries; i++)
      file_index->Fill();

    
    fin_tree->Write();
    fin_file->Close();
    std::cout<<ansi_green<<"[INDEXED] "<<ansi_reset<<Setup->vMCs.at(ifile)<<std::endl;
  }

  for(Int_t ifile = 0; ifile<(Int_t)Setup->vMCs.size(); ifile++)
    Setup->vMCs.at(ifile) = Form("IndexedFiles/Indexed_file_from_original_file_%i.root",ifile);

  return;

}



#endif
