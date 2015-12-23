///###########################################################################################################
///=========================[ Example of a main program using FastME module ]=================================
///==============================[ Code Author: Miqueias M. de Almeida ]======================================
///###########################################################################################################

#include <string>
#include <vector>

///FastME Module
#include "FastME.C"

void main_code(void){

  string path = "CMS_Ntuples/FastME_format/";
  string Data_Path(path+"4l_data.root");
  //string Data_Path(path+"4l_SM_Higgs126_1.root");
  //string Data_Path(path+"4l_ggZZ_plus_qqZZ.root"); 
 
  vector<string> MCs, MC_Names;

  //MCs.push_back(path+"4l_SM_Higgs126_1.root");
  MCs.push_back(path+"4l_SM_Higgs126_2.root");
  //MCs.push_back(path+"4l_ggZZ_1.root");
  MCs.push_back(path+"4l_ggZZ_2.root");
  //MCs.push_back(path+"4l_qqZZ_1.root");
  MCs.push_back(path+"4l_qqZZ_2.root");
  
  MC_Names.push_back("SM_Higgs126");
  MC_Names.push_back("ggZZ");
  MC_Names.push_back("qqZZ");


  ///Parameters (from left to right):
  ///1. Address to Data sample;
  ///2. Vector with address of MC samples;
  ///3.	Name of tree containing	the events;
  ///4. Vector with names to each MC type;
  ///5. Number of final state particles;
  ///6. Number of cores to be used;
  ///7. Name of the output file to store FastME analysis results.
  ///8. Verbose (0= quiet, 1= main info, 2= everything)

  FastME(Data_Path, MCs, "Higgs", MC_Names, 4, 3, "fme_results_tmp",2);

  ///FastME function will build a tree named "FastME" containing the discriminant value for each MC background 
  ///sample and the combined value when all background are considered. Also, it will store in the file 
  ///"fme_results_..." a TH2D histogram for the minimum distances found between the sample test and each MC 
  ///sample (signal & backgrounds).


}
