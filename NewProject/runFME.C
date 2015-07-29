///Interface to FastME
///Author: Miquéias M. Almeida

//Fast Matrix Element modules
#include "FastME.cxx"
#include "FastME.h"

#include <fstream>
#include <TString.h>
using namespace std;

int runFME(ifstream &Input){
  
  TString final_state, model, out_name, file_data_name, file_mcsig_name, file_mcbkg_name, tree_name, branch_name;
  TString info;

  Input >> info;
  do{
    if(info == "") continue;
    
    else if(info == "#OutName"){
      do{
	 Input >> info;
	 if(info == "#FinalState") break;
	 if(info.First("%")==0) continue;
	 out_name = info;
      }while(true);
    }
    else if(info == "#FinalState"){
      do{
	 Input >> info;
	 if(info == "#Model") break;
	 if(info.First("%")==0) continue;
	 final_state = info;
      }while(true);
    }
    else if (info == "#Model"){
      do{
	 Input >> info;
	 if(info == "#DataPath") break;
	 if(info.First("%")==0) continue;
	 model = info;
      }while(true);
    }
    else if(info == "#DataPath"){
      do{
	 Input >> info;
	 if(info == "#SigPath") break;
	 if(info.First("%")==0) continue;
	 file_data_name = info;
      }while(true);
    }
    else if(info == "#SigPath"){
      do{
	 Input >> info;
	 if(info == "#BkgPath") break;
	 if(info.First("%")==0) continue;
	 file_mcsig_name = info;
      }while(true);
    }
    else if(info == "#BkgPath"){
      do{
	 Input >> info;
	 if(info == "#TreeName") break;
	 if(info.First("%")==0) continue;
	 file_mcbkg_name = info;
      }while(true);
    }
    else if(info == "#TreeName"){
      do{
	 Input >> info;
	 if(info == "#BranchName") break;
	 if(info.First("%")==0) continue;
	 tree_name = info;
      }while(true);
    }
    else if(info == "#BranchName"){
      do{
	 Input >> info;
	 if(info == "#fim") break;
	 if(info.First("%")==0) continue;
	 branch_name = info;
      }while(true);
    }
  }while(info != "#fim");
  
  ///Instatiating and activing FastME analysis
  FME *mFME = new FME();
  mFME->launchFME(final_state,model,out_name,file_data_name,file_mcsig_name,file_mcbkg_name,tree_name,branch_name);

  return 0;
}