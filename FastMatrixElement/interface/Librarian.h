#ifndef Librarian_h
#define Librarian_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

void Indexer(FmeSetup Setup){

  ///Inserts a branch into the file to store its index in the config file
  std::cout<<":: Indexing files: "<<std::endl;  

  for(Int_t ifile = 0; ifile < (Int_t)Setup.vMCs.size(); ifile++){

    TFile *org_file = new TFile( (TString)Setup.vMCs.at(ifile), "update" );
    TTree *org_tree = (TTree*)org_file->Get(Setup.TTreeName);
    Int_t McFileIndex = ifile;
    TBranch *file_index = org_tree->Branch("McFileIndex",&McFileIndex,"McFileIndex/I");
    Int_t Nentries = org_tree->GetEntries();

    for(Int_t i = 0; i < Nentries; i++)
      file_index->Fill();
    org_tree->Write();

    std::cout<<ansi_green<<"[READY] "<<ansi_reset<<org_file->GetName()<<std::endl;
    delete org_file;
  }


  return;

}


#endif
