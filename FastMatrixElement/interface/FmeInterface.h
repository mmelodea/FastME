///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Fast Matrix Element Interface Manager ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/FileFormater.h"
//#include "FastMatrixElement/FastMatrixElement/interface/ShowParticles.h"


#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <stdlib.h>



///Contains the menu of available commands
static std::string help = "-help", nc = "-c", ff = "-f", fa = "-a", sp = "-s";
void Helper(void){
  std::cout<<"Usage: fastme [commands] config_file"<<std::endl;
  std::cout<<"Commands:"<<std::endl;
  std::cout<<"\t-c\t\tInform how many cores are available in the machine"<<std::endl;
  std::cout<<"\t-f\t\tConvert a general root file to FastME root file format"<<std::endl;
  std::cout<<"\t-a\t\tMake the FastME analysis over events"<<std::endl;
  std::cout<<"\t-s\t\tDisplay the particles disposition on FastME phase space"<<std::endl;  
  std::cout<<"More info access github..."<<std::endl;

  return;
}


///Return how many cores are available in the machine
void FindCores(){
  std::cout<<Form("cores available: %i",system("nproc"))<<endl;
  
  return;
}


///Read input file and convert to program format
void ConfigReader(std::string UserConfig, FmeSetup *Setup){    
  ///Define variables used in the analysis
  TString Data_Path;
  std::vector<TString> MC_Names;
  
  ///______________ Extract the need info from txt file _______________________________________________________________
  ///Opens the config file to get the user configuration
  std::ifstream inFile(UserConfig.c_str());
  if(!inFile){
    std::cout<<"Error in file!";
    throw std::exception();
  }
  int nkeys=0;
  std::string line;
  while( getline(inFile,line) ){
    if(line.find("#") != std::string::npos) continue;
    for(int k=0; k<(int)fme_keywords.size(); k++){
      if(line.find(fme_keywords[k]) != std::string::npos){
	nkeys++;
        line.erase(line.begin(),line.begin()+ksize[k]);
        if(k==0 || k==1 || k==2 || k==3 || k==8 || k==9 || k==13 || k==14 || k==15 || k==16 || k==17 || k==19)
	  std::cout <<":: "<< fme_keywords[k] <<"\t\t\t\t"<< line << std::endl;
	else std::cout <<":: "<< fme_keywords[k] <<"\t\t\t"<< line << std::endl;
        if(fme_keywords[k] == 		"data_path") 	Data_Path = line;
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
	if(fme_keywords[k] == 	   "n_fs_particles") 	Setup->NFSParticles = stoi(line);
	if(fme_keywords[k] == 	"flavor_constraint")	Setup->SetFlavorConstraint = line;
	if(fme_keywords[k] == 		  "n_cores")	Setup->NCores = stoi(line);
	if(fme_keywords[k] == 	       "data_limit")	Setup->DTLimit = stoi(line);
	if(fme_keywords[k] == 		 "mc_limit")	Setup->MCLimit = stof(line);
        if(fme_keywords[k] == 		"scale_dPt")	Setup->ScaledPt = stof(line);
        if(fme_keywords[k] == 	       "scale_dEta")	Setup->ScaledEta = stof(line);
	if(fme_keywords[k] == 	     "scale_method")	Setup->ScaleMethod = line;
	if(fme_keywords[k] == 	    "verbose_level")	Setup->Verbose = stoi(line);
      }
    }
  }
  if(nkeys < ((int)fme_keywords.size()-2)){
    std::cout<<"Missing key-word! Check your input file!";
    throw std::exception();
  }
  ///__________________________________________________________________________________________________________________
  
  ///Getting some numbers
  TFile *fData = TFile::Open(Data_Path);
  TTreeReader tmpReader1(Setup->TTreeName,fData);
  Int_t nData = tmpReader1.GetEntries(true);
  if(Setup->DTLimit != -1 && Setup->DTLimit <= nData)
    nData = Setup->DTLimit;

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
  std::cout<<"______________________________________________________________________________________________"<<std::endl;
  std::cout<< Form(":: Data Events:      %i",nData) <<std::endl;
  std::cout<< Form(":: MC Samples:       %i\t[",N_MC); 
            for(Int_t ne=0; ne<N_MC; ne++){
              if( ne<(N_MC-1) ) std::cout<< ne << "= " << NMCEV[ne] << ",  ";
              if( ne==(N_MC-1) ) std::cout<< ne << "= " << NMCEV[ne] << "]" <<std::endl;
            }
  std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;
  ///--------------------------------------------------------------------------------------------------------------
  
  Setup->DataFile = fData;
  Setup->NData = nData;
  Setup->NMCT = N_MCT;
  
  return;
}



///--------------------- Interface manager ----------------------
int FmeInterface(char *argv[], FmeSetup *USetup){
       if(argv[1] == help)	Helper();
  else if(argv[1] == nc)	FindCores();
  else if(argv[1] == ff)	FileFormater((std::string)argv[2]);
  else if(argv[1] == fa)	ConfigReader((std::string)argv[2], USetup);
  //else if(argv[1] == sp)	ShowParticles();
  else{
    std::cout<<"[ERROR] Invalid command '"<<argv[1]<<"'"<<std::endl;
    std::cout<<"These are the available commands:"<<std::endl;
    Helper();
    return -1;
  }
  
  return 0;
}
