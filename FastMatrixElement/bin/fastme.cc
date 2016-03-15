///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::: FastME executable ::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/FmeInterface.h"
#include "FastMatrixElement/FastMatrixElement/interface/ComputePhsDR.h"
#include "FastMatrixElement/FastMatrixElement/interface/Discriminant.h"


#include <iostream>
#include <string>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TSystem.h>


int main(int argc, char *argv[]){

  ///Instantiate the needed variables
  FmeSetup setup;
  TTree *rtree, *ftree;

  ///-------------------------------------------------------------------------------------------------------------------
  std::cout<<"\n\n";
  std::cout<<"=============================================================================================="<<std::endl;
  std::cout<<"::::::::::::::::::::::::::[ Fast Matrix Element Analysis Started ]::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<"=============================================================================================="<<std::endl;
  std::cout<<":: [Your input file: "<< argv[2] <<"]"<<std::endl;
  
  ///Calls interface to get analysis setup
  FmeInterface(argv, &setup);

  ///User can abort analysis
  std::cout<<" Input file correct?(y/n) ";
  std::string aws3;
  std::cin >> aws3;
  if(aws3 == "n") return -1;  
  
  ///Calls PhsDrComputer to compute events distance
  rtree = ComputePhsDR(setup);
  
  ///Calls Discriminator
  ftree = Discriminant(rtree, setup);

  gSystem->Exec("mkdir -p "+setup.OutPath);
  TString resulting_file = setup.OutPath+"/"+setup.OutName+".root";
  TFile *ffile = new TFile(resulting_file,"recreate");
  ftree->Write();
  ffile->Close();
  
  ///-------------------------------------------------------------------------------------------------------------------
  std::cout<<"=============================================================================================="<<std::endl;
  std::cout<<":::::::::::::::::::::::::[ Fast Matrix Element Analysis Finalized ]:::::::::::::::::::::::::::"<<std::endl;
  //std::cout<<":: [Analysis file saved: "<<resulting_file<<"]"<<std::endl;
  std::cout<<"=============================================================================================="<<std::endl;
  std::cout<<"\n\n";

  return 0;
}
