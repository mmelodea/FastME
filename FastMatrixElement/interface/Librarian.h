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
  std::cout<<"\n>>> Actioning Librarian to handle the input files <<<"<<std::endl;

  ///Create a folder to keep data to be used by FastME
  ///Only the important branches and the defined events quantity are keept
  gSystem->Exec("mkdir FME_USAGE");
  gSystem->Exec("mkdir FME_USAGE/DATA");
  gSystem->Exec("mkdir FME_USAGE/MC_TEMPLATES");


  ///Scale the number of events by Cross Section
  Int_t max_sig_xs_data_ev = 0, max_bkg_xs_data_ev = 0;
  float max_sig_xs_data = 0, max_bkg_xs_data = 0;
  Int_t max_sig_xs_mc_ev = 0, max_bkg_xs_mc_ev = 0;
  float max_sig_xs_mc = 0, max_bkg_xs_mc = 0;
  if(Setup->XSScale == "true"){    
    std::vector<int> SigData = Setup->SigData;
    std::vector<int> BkgData = Setup->BkgData;
    std::vector<float> DataXS = Setup->DataXS;
    std::vector<float> McXS = Setup->McXS;
    
    if(DataXS.size() != 0){
      std::cout<<"Computing scaling factors for each data file..."<<std::endl;
      int max_sig_xs_id = -1, max_bkg_xs_id = -1;
      for(Int_t ixs=0; ixs<(Int_t)DataXS.size(); ixs++){
	if(ixs < SigData.size() && DataXS[ixs] > max_sig_xs_data){
	  max_sig_xs_data = DataXS[ixs];
	  max_sig_xs_id = ixs;
	}
	if(ixs >= SigData.size() && DataXS[ixs] > max_bkg_xs_data){
	  max_bkg_xs_data = DataXS[ixs];
	  max_bkg_xs_id = ixs;
	}
      }
      
      TFile *tmpfile = TFile::Open(Setup->vDatas[max_sig_xs_id]);
      TTree *tmptree = (TTree*)tmpfile->Get(Setup->TTreeName);
      max_sig_xs_data_ev = tmptree->GetEntries();
      tmpfile->Close();

      TFile *tmpfile = TFile::Open(Setup->vDatas[max_bkg_xs_id]);
      TTree *tmptree = (TTree*)tmpfile->Get(Setup->TTreeName);
      max_bkg_xs_data_ev = tmptree->GetEntries();
      tmpfile->Close();      
    }

    if(McXS.size() != 0){
      std::cout<<"Computing scaling factors for each mc file..."<<std::endl;
      int max_sig_xs_id = -1, max_bkg_xs_id = -1;
      for(Int_t ixs=0; ixs<(Int_t)McXS.size(); ixs++){
	if(ixs < SigMC.size() && McXS[ixs] > max_sig_xs_mc){
	  max_sig_xs_mc = McXS[ixs];
	  max_sig_xs_id = ixs;
	}
	if(ixs >= SigMC.size() && McXS[ixs] > max_bkg_xs_mc){
	  max_bkg_xs_mc = McXS[ixs];
	  max_bkg_xs_id = ixs;
	}
      }
      
      TFile *tmpfile = TFile::Open(Setup->vMCs[max_sig_xs_id]);
      TTree *tmptree = (TTree*)tmpfile->Get(Setup->TTreeName);
      max_sig_xs_mc_ev = tmptree->GetEntries();
      tmpfile->Close();

      TFile *tmpfile = TFile::Open(Setup->vDatas[max_bkg_xs_id]);
      TTree *tmptree = (TTree*)tmpfile->Get(Setup->TTreeName);
      max_bkg_xs_mc_ev = tmptree->GetEntries();
      tmpfile->Close();      
    }
    
  }
  
  ///Preparing the Data files
  std::cout<<":: Reducing data files: "<<std::endl;  
  std::cout<<"::------------------------------------------------------------::"<<std::endl;
  Int_t nDataFiles = Setup->vDatas.size();
  for(Int_t ifile = 0; ifile < nDataFiles; ++ifile){
    Int_t limit;
    if(ifile < Setup->SigData.size()) limit = (Setup->DataXS[ifile]/max_sig_xs_data)*max_sig_xs_data_ev;
    else limit = (Setup->DataXS[ifile]/max_bkg_xs_data)*max_bkg_xs_data_ev;

    TFile *org_file = TFile::Open( (TString)Setup->vDatas[ifile] );
    TTree *org_tree = (TTree*)org_file->Get(Setup->TTreeName);
    org_tree->SetBranchStatus("*",kFALSE);
    org_tree->SetBranchStatus(Setup->IdBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->PtBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->EtaBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->PhiBranch,kTRUE);
    if(limit > org_tree->GetEntries()){
      std::cout<<"File doesn't have envents enough to scale! Using all them."<<std::endl;
      limit = org_tree->GetEntries();
    }

    TFile *fin_file = new TFile( Form("FME_USAGE/DATA/Reduced_file_from_original_data_file_%i.root",ifile), "recreate" );
    TTree *fin_tree = org_tree->CloneTree(limit);
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
    Int_t limit;
    if(ifile < Setup->SigMC.size()) limit = (Setup->McXS[ifile]/max_sig_xs_mc)*max_sig_xs_mc_ev;
    else limit = (Setup->McXS[ifile]/max_bkg_xs_mc)*max_bkg_xs_mc_ev;

    TFile *org_file = TFile::Open( (TString)Setup->vMCs[ifile] );
    TTree *org_tree = (TTree*)org_file->Get(Setup->TTreeName);
    org_tree->SetBranchStatus("*",kFALSE);
    org_tree->SetBranchStatus(Setup->IdBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->PtBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->EtaBranch,kTRUE);
    org_tree->SetBranchStatus(Setup->PhiBranch,kTRUE);
    if(limit > org_tree->GetEntries()){
      std::cout<<"File doesn't have envents enough to scale! Using all them."<<std::endl;
      limit = org_tree->GetEntries();
    }

    TFile *fin_file = new TFile( Form("FME_USAGE/MC_TEMPLATES/Indexed_file_from_original_MC_file_%i.root",ifile), "recreate" );
    TTree *fin_tree = org_tree->CloneTree(limit);
    org_file->Close();

    Int_t McFileIndex = ifile;
    TBranch *file_index = fin_tree->Branch("McFileIndex",&McFileIndex,"McFileIndex/I");
    for(Int_t i = 0; i < limit; ++i) file_index->Fill();

    fin_tree->Write();
    fin_file->Close();
    std::cout<<ansi_green<<"[INDEXED] "<<ansi_reset<<Setup->vMCs[ifile]<<std::endl;
    Setup->vMCs[ifile] = Form("FME_USAGE/MC_TEMPLATES/Indexed_file_from_original_MC_file_%i.root",ifile);
  }




  return;

}



#endif
