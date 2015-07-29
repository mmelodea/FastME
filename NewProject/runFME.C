///Interface to FastME
///Author: Miquéias M. Almeida

//Fast Matrix Element modules
#include "FastME.cxx"
#include "FastME.h"

#include <fstream>
#include <TString.h>
using namespace std;

int runFME(ifstream &Input){
  
  TString model, out_name, file_data_name, file_mcsig_name, file_mcbkg_name, tree_name, branch_name;
  TString info;
  do{
    Input >> info;
    if(info == "") continue;
    if(info == "#Data_path"){
      Input >> info;
      file_data_name = info;
    }
    else if(info == "#Sig_path"){
      Input >> info;
      file_mcsig_name = info;
    }
    else if(info == "#Bkg_path"){
      Input >> info;
      file_mcbkg_name = info;
    }
    else if(info == "#Tree_name"){
      Input >> info;
      tree_name = info;
    }
    else if(info == "#Branch_name"){
      Input >> info;
      branch_name = info;
    }
    else if(info == "#Out_name"){
      Input >> info;
      out_name = info;
    }
    else if(info == "#Model"){
      Input >> info;
      model = info;
    }
  }while(info != "#fim");
  
  ///Instatiating and activing FastME analysis
  FME *mFME = new FME();
  mFME->launchFME(model,out_name,file_data_name,file_mcsig_name,file_mcbkg_name,tree_name,branch_name);

  return 0;
}