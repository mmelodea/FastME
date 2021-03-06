///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Fast Matrix Element Interface Manager ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Cortana_h
#define Cortana_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/InitScreen.h"
#include "FastMatrixElement/FastMatrixElement/interface/Librarian.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <iomanip>




///Test file existence
void file_check(std::string file_address){
  std::ifstream exists(file_address.c_str());
  if(!exists){
    std::cout<<"File "<<file_address<<" not found!"<<std::endl;
    throw std::exception();
  }

  return;
}



///Contains the menu of available commands
static std::string help = "-help", sl = "-quiet", normal = "noprint", nc = "-c", fa = "-a", pr = "-d", ge = "-g";
void Helper(void){
  std::cout<<"Usage: fastme [commands] config_file [flux options]"<<std::endl;
  std::cout<<"Commands:"<<std::endl;
  std::cout<<"\t-c\t\tInform how many cores are available in the machine"<<std::endl;
  std::cout<<"\t-a\t\tMake the FastME analysis over events"<<std::endl;
  std::cout<<"\t-d\t\tCompute discriminant and produce plots from the FastME results"<<std::endl;
  std::cout<<"\t-g\t\tGenerate events using the minimum distance method"<<std::endl;
  std::cout<<"Flux options:"<<std::endl;
  std::cout<<"\t-quiet\t\tRun program without show the config file and avoid user check config file"<<std::endl;
  std::cout<<"\nFor more info access https://github.com/mmelodea/FastMatrixElement"<<std::endl;

  return;
}



///Return how many cores are available in the machine
void FindCores(){
  system("nproc");
  std::cout<<"Cores available"<<std::endl;
  std::cout<<ansi_yellow<<"[Warning] Be carefull, do not set a number of cores higher than available!"<<ansi_reset<<std::endl;
  
  return;
}


///Read input file and convert to program format
void ConfigReader(std::string UserConfig, FmeSetup *Setup, std::string command, std::string run_mode = normal){

  ///Print FastME name over screen
  if(run_mode != sl){
    InitScreen();
  }

  ///Define variables used in the analysis
  TString Data_Path;
  
  ///______________ Extract the need info from txt file _______________________________________________________________
  ///Opens the config file to get the user configuration
  std::ifstream inFile(UserConfig.c_str());
  if(!inFile){
    std::cout<<ansi_red<<"[ERROR]"<<ansi_reset<<" File could not be openned!";
    throw std::exception();
  }

  std::cout<<"\n:: ["<<ansi_yellow<<"Your input file..."<<ansi_reset<<"]"<<std::endl;
  int nkeys=0;
  std::string line;
  while( getline(inFile,line) ){
    if(line.find("#") != std::string::npos) continue;
    for(int k=0; k<(int)fme_keywords.size(); k++){
      if(line.find(fme_keywords.at(k)) != std::string::npos){
	nkeys++;
        line.erase(line.begin(),line.begin()+ksize.at(k));
	if(run_mode == normal){
	  if(command == pr){
            if(fme_keywords.at(k) == "fme_file")
	      std::cout << ":: " << ansi_cyan << fme_keywords.at(k) << ":  " << ansi_reset << line << std::endl;
          }
          else{
    	      std::cout << ":: " << ansi_cyan << fme_keywords.at(k) << ":  " << ansi_reset << line << std::endl;
          }
	}
        if(fme_keywords.at(k) ==             "data_path")	Setup->vDatas.push_back(line);
	if(fme_keywords.at(k) ==               "mc_path")	Setup->vMCs.push_back(line);
	if(fme_keywords.at(k) ==             "tree_name")	Setup->TTreeName = line;
	if(fme_keywords.at(k) ==        "id_branch_name")	Setup->IdBranch = line;
	if(fme_keywords.at(k) ==        "pt_branch_name")	Setup->PtBranch = line;
	if(fme_keywords.at(k) ==       "eta_branch_name")	Setup->EtaBranch = line;
	if(fme_keywords.at(k) ==       "phi_branch_name")	Setup->PhiBranch = line;
	if(fme_keywords.at(k) ==          "outfile_name")	Setup->OutName = line;
	if(fme_keywords.at(k) ==          "outfile_path")	Setup->OutPath = line;
	if(fme_keywords.at(k) ==         "phs_dr_method")	Setup->PhSDrMethod = line;
	if(fme_keywords.at(k) ==     "flavor_constraint")	Setup->SetFlavorConstraint = line;
	if(fme_keywords.at(k) ==               "n_cores")	Setup->NCores = stoi(line);
	if(fme_keywords.at(k) ==            "data_limit")	Setup->DTLimit = stoi(line);
	if(fme_keywords.at(k) ==              "mc_limit")	Setup->MCLimit = stof(line);
        if(fme_keywords.at(k) ==             "scale_dPt")	Setup->ScaledPt = stof(line);
        if(fme_keywords.at(k) ==            "scale_dEta")	Setup->ScaledEta = stof(line);
        if(fme_keywords.at(k) ==            "scale_dPhi")	Setup->ScaledPhi = stof(line);
	if(fme_keywords.at(k) ==              "fme_file")	Setup->FmeFile = line;
        if(fme_keywords.at(k) ==              "sig_data")       Setup->SigData.push_back(stoi(line));
        if(fme_keywords.at(k) ==              "bkg_data")       Setup->BkgData.push_back(stoi(line));
        if(fme_keywords.at(k) ==                "sig_mc")       Setup->SigMC.push_back(stoi(line));
        if(fme_keywords.at(k) ==                "bkg_mc")       Setup->BkgMC.push_back(stoi(line));
        if(fme_keywords.at(k) ==              "xs_scale")       Setup->XSScale = line;
        if(fme_keywords.at(k) ==               "data_xs")       Setup->DataXS.push_back(stof(line));
        if(fme_keywords.at(k) ==                 "mc_xs")       Setup->McXS.push_back(stof(line));
	if(fme_keywords.at(k) ==            "gen_factor")	Setup->GenFactor = stoi(line);
	if(fme_keywords.at(k) ==               "dist_pt")	Setup->DistPt = line;
	if(fme_keywords.at(k) ==              "dist_eta")	Setup->DistEta = line;
	if(fme_keywords.at(k) ==              "dist_phi")	Setup->DistPhi = line;
	if(fme_keywords.at(k) ==         "gaussian_mean")	Setup->GaussianMean = stof(line);
        if(fme_keywords.at(k) ==         "sdr_condition")       Setup->SDrCondition = stof(line);
        if(fme_keywords.at(k) ==         "bdr_condition")       Setup->BDrCondition = stof(line);
        if(fme_keywords.at(k) ==        "max_gen_trials")       Setup->MaxGenTrials = stoi(line);
	if(fme_keywords.at(k) ==          "out_gen_name")       Setup->OutGenName = line;
	if(fme_keywords.at(k) ==         "verbose_level")	Setup->Verbose = stoi(line);
      }
    }
  }//While's end
  //if(nkeys < ((int)fme_keywords.size()-2)){
    //std::cout<<ansi_red<<"[ERROR]"<<ansi_reset<<" Missing key-word! Check your input file!";
    //throw std::exception();
  //}
  ///__________________________________________________________________________________________________________________

  

  if(command == ge) return;
  
  ///Getting some info
  std::cout<<"\n:: "<<ansi_yellow<<"Checking inputs..."<<ansi_reset<<std::endl;
  

  const Int_t N_DT = Setup->vDatas.size();
  Int_t NDATA[N_DT];
  Int_t nDtEv = 0;

  if(command != pr)
    for(Int_t nd=0; nd<(Int_t)Setup->vDatas.size(); nd++){
      file_check(Setup->vDatas.at(nd));
      TFile *fData = TFile::Open((TString)Setup->vDatas.at(nd));
      TTreeReader tmpReader1(Setup->TTreeName,fData);
      Int_t nData = tmpReader1.GetEntries(true);
      if(Setup->DTLimit != -1 && Setup->DTLimit <= nData)
      nData = Setup->DTLimit;
      NDATA[nd] = nData;
      nDtEv += nData;
      fData->Close();
    }


  const Int_t N_MC = Setup->vMCs.size();
  Int_t NMCEV[N_MC];

  if(command != pr)
    for(Int_t ne=0; ne<N_MC; ne++){
      file_check(Setup->vMCs.at(ne));
      TFile *fmc = TFile::Open((TString)Setup->vMCs.at(ne));
      TTreeReader tmpReader2(Setup->TTreeName,fmc);
      if(Setup->MCLimit == -1)
        NMCEV[ne] = tmpReader2.GetEntries(true);
      if(Setup->MCLimit != -1 && Setup->MCLimit < 1)
        NMCEV[ne] = (Int_t)(Setup->MCLimit*tmpReader2.GetEntries(true));
      if(Setup->MCLimit != -1 && Setup->MCLimit >= 1)
        NMCEV[ne] = Setup->MCLimit;
      fmc->Close();
    }
  if(N_MC == 0){
    std::cout<<ansi_red<<"[ERROR]"<<ansi_reset<<" None MC template detected! Please, revise your input file!";
    throw std::exception();
  }

  if(run_mode == normal && command != pr){
    Int_t n_dataev = 0, n_mcev = 0;
    std::cout<<ansi_yellow<<"______________________________________________________________________________________________"<<ansi_reset<<std::endl;
    std::cout<< Form(":: Data Samples:     %i\t[",N_DT);
    for(Int_t nd = 0; nd < N_DT; ++nd){
      if( nd<(N_DT-1) ) std::cout<< nd << "= " << NDATA[nd] << ",  ";
      if( nd==(N_DT-1) ) std::cout<< nd << "= " << NDATA[nd] << "]" <<std::endl;
      n_dataev += NDATA[nd];
    }
    std::cout<< Form(":: MC Samples:       %i\t[",N_MC); 
    for(Int_t ne = 0; ne < N_MC; ++ne){
      if( ne<(N_MC-1) ) std::cout<< ne << "= " << NMCEV[ne] << ",  ";
      if( ne==(N_MC-1) ) std::cout<< ne << "= " << NMCEV[ne] << "]" <<std::endl;
      n_mcev += NMCEV[ne];
    }
    std::cout<< Form(":: Going to analyze %i through %i Monte Carlo events", n_dataev, n_mcev) <<std::endl;
    std::cout<<ansi_yellow<<"----------------------------------------------------------------------------------------------"<<ansi_reset<<std::endl;
  }
  ///--------------------------------------------------------------------------------------------------------------

  Setup->NData = nDtEv;

  ///Get files order and insert a branch inside them to handle in different CPU cores
  ///Indexer function defined into Librarian header
  if(command != pr)  Indexer(Setup);

  
  return;
}


#endif
