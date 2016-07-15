///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Fast Matrix Element Interface Manager ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef FmeInterface_h
#define FmeInterface_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/InitScreen.h"

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



///Contains the menu of available commands
static std::string help = "-help", sl = "-quiet", normal = "noprint", nc = "-c", fa = "-a", pr = "-p", sp = "-s";
void Helper(void){
  std::cout<<"Usage: fastme [commands] config_file [flux options]"<<std::endl;
  std::cout<<"Commands:"<<std::endl;
  std::cout<<"\t-c\t\tInform how many cores are available in the machine"<<std::endl;
  std::cout<<"\t-a\t\tMake the FastME analysis over events"<<std::endl;
  std::cout<<"\t-p\t\tMake plots from the FastME results"<<std::endl;
  std::cout<<"\t-s\t\tDisplay the particles disposition on FastME phase space"<<ansi_yellow<<" (to be implemented)"<<ansi_reset<<std::endl;  
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
  std::vector<TString> MC_Names;
  
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
      if(line.find(fme_keywords[k]) != std::string::npos){
	nkeys++;
        line.erase(line.begin(),line.begin()+ksize[k]);
	if(run_mode == normal){
	  if(command == pr){
            if(fme_keywords[k] == "fme_files")
	      std::cout << ":: " << ansi_cyan << fme_keywords[k] << ":  " << ansi_reset << line << std::endl;
          }
          else{
    	      std::cout << ":: " << ansi_cyan << fme_keywords[k] << ":  " << ansi_reset << line << std::endl;
          }
	}
        if(fme_keywords[k] == 		"data_path") 	Setup->vDatas.push_back(line);
	if(fme_keywords[k] == 		  "mc_path")	Setup->vMCs.push_back(line);
	if(fme_keywords[k] == 		  "mc_name")	MC_Names.push_back(line);
	if(fme_keywords[k] == 		"tree_name")	Setup->TTreeName = line;
	if(fme_keywords[k] == "mc_type_branch_name") 	Setup->McTypeBranch = line;
	if(fme_keywords[k] == 	   "id_branch_name") 	Setup->IdBranch = line;
	if(fme_keywords[k] == 	   "pt_branch_name") 	Setup->PtBranch = line;
	if(fme_keywords[k] == 	  "eta_branch_name") 	Setup->EtaBranch = line;
	if(fme_keywords[k] == 	     "outfile_name")	Setup->OutName = line;
	if(fme_keywords[k] == 	     "outfile_path")	Setup->OutPath = line;
	if(fme_keywords[k] == 	    "phs_dr_method")	Setup->PhSDrMethod = line;
	if(fme_keywords[k] == 	"flavor_constraint")	Setup->SetFlavorConstraint = line;
	if(fme_keywords[k] == 		  "n_cores")	Setup->NCores = stoi(line);
	if(fme_keywords[k] == 	       "data_limit")	Setup->DTLimit = stoi(line);
	if(fme_keywords[k] == 		 "mc_limit")	Setup->MCLimit = stof(line);
        if(fme_keywords[k] == 		"scale_dPt")	Setup->ScaledPt = stof(line);
        if(fme_keywords[k] == 	       "scale_dEta")	Setup->ScaledEta = stof(line);
	if(fme_keywords[k] == 	     	"fme_files")	Setup->FmeFiles.push_back(line);
        if(fme_keywords[k] ==        "storePhSTree")    Setup->StorePhSTree = line;
	if(fme_keywords[k] == 	    "verbose_level")	Setup->Verbose = stoi(line);
      }
    }
  }
  if(nkeys < ((int)fme_keywords.size()-2)){
    std::cout<<ansi_red<<"[ERROR]"<<ansi_reset<<" Missing key-word! Check your input file!";
    throw std::exception();
  }
  ///__________________________________________________________________________________________________________________
  
  ///Getting some numbers
  std::cout<<"\n:: "<<ansi_yellow<<"Checking inputs..."<<ansi_reset<<std::endl;
  const Int_t N_DT = Setup->vDatas.size();
  Int_t NDATA[N_DT];
  for(Int_t nd=0; nd<(Int_t)Setup->vDatas.size(); nd++){
    TFile *fData = TFile::Open((TString)Setup->vDatas[nd]);
    TTreeReader tmpReader1(Setup->TTreeName,fData);
    Int_t nData = tmpReader1.GetEntries(true);
    if(Setup->DTLimit != -1 && Setup->DTLimit <= nData)
    nData = Setup->DTLimit;

    NDATA[nd] = nData;
    fData->Close();
  }

  const Int_t N_MCT = MC_Names.size();
  const Int_t N_MC = Setup->vMCs.size();
  Int_t NMCEV[N_MC];
  for(Int_t ne=0; ne<N_MC; ne++){
    TFile *fmc = TFile::Open((TString)Setup->vMCs[ne]);
    TTreeReader tmpReader2(Setup->TTreeName,fmc);
    if(Setup->MCLimit == -1)
      NMCEV[ne] = tmpReader2.GetEntries(true);
    if(Setup->MCLimit != -1 && Setup->MCLimit < 1)
      NMCEV[ne] = (Int_t)(Setup->MCLimit*tmpReader2.GetEntries(true));
    if(Setup->MCLimit != -1 && Setup->MCLimit >= 1)
      NMCEV[ne] = Setup->MCLimit;
    fmc->Close();
  }
  if(run_mode == normal && command != pr){
    std::cout<<ansi_yellow<<"______________________________________________________________________________________________"<<ansi_reset<<std::endl;
    std::cout<< Form(":: Data Samples:     %i\t[",N_DT);
    for(Int_t nd=0; nd<N_DT; nd++){
      if( nd<(N_DT-1) ) std::cout<< nd << "= " << NDATA[nd] << ",  ";
      if( nd==(N_DT-1) ) std::cout<< nd << "= " << NDATA[nd] << "]" <<std::endl;
    }
    std::cout<< Form(":: MC Samples:       %i\t[",N_MC); 
    for(Int_t ne=0; ne<N_MC; ne++){
      if( ne<(N_MC-1) ) std::cout<< ne << "= " << NMCEV[ne] << ",  ";
      if( ne==(N_MC-1) ) std::cout<< ne << "= " << NMCEV[ne] << "]" <<std::endl;
    }
    std::cout<<ansi_yellow<<"----------------------------------------------------------------------------------------------"<<ansi_reset<<std::endl;
  }
  ///--------------------------------------------------------------------------------------------------------------
  
  
  Setup->NData = N_DT;
  Setup->NMCT = N_MCT;
  
  return;
}


#endif
