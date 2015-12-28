///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::[ Example on How to Use FastME Module ]:::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::[ Code Author: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#include <string>
#include <vector>
#include "FastME.C"

void main_code(void){
  string path = "../CMS_Ntuples/FastME_format/";
  //string Data_Path(path+"4l_data.root");
  //string Data_Path(path+"4l_SM_Higgs126_1.root");
  string Data_Path(path+"4l_ggZZ_plus_qqZZ.root"); 
 
  vector<string> MCs;

  //MCs.push_back(path+"4l_SM_Higgs126_1.root");
  MCs.push_back(path+"4l_SM_Higgs126_2.root");
  //MCs.push_back(path+"4l_ggZZ_1.root");
  MCs.push_back(path+"4l_ggZZ_2.root");
  //MCs.push_back(path+"4l_qqZZ_1.root");
  MCs.push_back(path+"4l_qqZZ_2.root");

  ///Parameters:
  ///1. Address to Data sample;
  ///2. Vector with address of MC samples;
  //FastME(Data_Path, MCs);
  FastME();

  ///Here just these parameters are given. However there's a file called "fme_config.dat" at FME_module, in which 
  ///the user can specify many other configs (like samples info - tree name, branches - comparison method, MC names, name of output, etc). This is done
  ///to avoid delay: compile code and pass samples by argument slow down the analysis in half of the time spent if
  ///the samples are uploaded to RAM memory. 
  ///Results are saved by default at FME_results
}
