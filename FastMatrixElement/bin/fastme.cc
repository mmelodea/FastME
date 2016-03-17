///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::[ Fast Matrix Element Executable ]::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/FmeInterface.h"
#include "FastMatrixElement/FastMatrixElement/interface/ComputePhsDR.h"
#include "FastMatrixElement/FastMatrixElement/interface/Discriminant.h"
#include "FastMatrixElement/FastMatrixElement/interface/StudyResults.h"
#include "FastMatrixElement/FastMatrixElement/interface/FileFormater.h"
//#include "FastMatrixElement/FastMatrixElement/interface/ShowParticles.h"


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

  ///Takes the config input file and converts it to FastME readable format
  if(argv[1] == fa)	ConfigReader((std::string)argv[2], &setup);

  ///--------------------- Interface manager ------------------------------
  ///To show usage
  if(argv[1] == help)	Helper();
  
  ///To show number of cores available
  else if(argv[1] == nc)	FindCores();
  
  ///To format an ntuple in different format
  else if(argv[1] == ff)	FileFormater(setup);
  
  ///To make the Fast Matrix Element analysis and compute discriminant
  else if(argv[1] == fa){
    std::cout<<"\n\n";
    std::cout<<"=============================================================================================="<<std::endl;
    std::cout<<"::::::::::::::::::::::::::[ Fast Matrix Element Analysis Started ]::::::::::::::::::::::::::::"<<std::endl;
    std::cout<<"=============================================================================================="<<std::endl;
    std::cout<<":: [Your input file: "<< argv[2] <<"]"<<std::endl;
    ///User can abort analysis if something is wrong
    std::cout<<" Input file correct?(y/n) ";
    std::string aws3;
    std::cin >> aws3;
    if(aws3 == "n") return -1;  

    ///Calls PhsDrComputer to compute events distance
    rtree = ComputePhsDR(setup);

    ///Calls Discriminator
    ftree = Discriminant(rtree, setup);

    ///Finalize results
    gSystem->Exec("mkdir -p "+setup.OutPath);
    TString resulting_file = setup.OutPath+"/"+setup.OutName+".root";
    TFile *ffile = new TFile(resulting_file,"recreate");
    ftree->Write();
    ffile->Close();
  
    ///-------------------------------------------------------------------------------------------------------------------
    std::cout<<"=============================================================================================="<<std::endl;
    std::cout<<":::::::::::::::::::::::::[ Fast Matrix Element Analysis Finalized ]:::::::::::::::::::::::::::"<<std::endl;
    std::cout<<":: [Analysis file saved: "<<resulting_file<<"]"<<std::endl;
    std::cout<<"=============================================================================================="<<std::endl;
    std::cout<<"\n\n";
  }

  ///Calls FastME analyzer (plot discriminants, ROC curve, events/discriminant value)
  else if(argv[1] == pr )	StudyResults(setup);
  
  ///Calls event display
  //else if(argv[1] == sp)	ShowParticles();
  
  ///Wrong commands gets error
  else{
    std::cout<<"[ERROR] Invalid command '"<<argv[1]<<"'"<<std::endl;
    std::cout<<"These are the available commands:"<<std::endl;
    Helper();
    return -1;
  }
  ///----------------------------------------------------------------------  
 
  return 0;
}
