///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::[ Fast Matrix Element Executable ]::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



///My Own headers
#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"
#include "FastMatrixElement/FastMatrixElement/interface/FmeInterface.h"
#include "FastMatrixElement/FastMatrixElement/interface/ComputePhsDR.h"
#include "FastMatrixElement/FastMatrixElement/interface/Discriminant.h"
#include "FastMatrixElement/FastMatrixElement/interface/StudyResults.h"
#include "FastMatrixElement/FastMatrixElement/interface/FileFormater.h"
#include "FastMatrixElement/FastMatrixElement/interface/ShowParticles.h"

//c++ headers
#include <iostream>
#include <string>

//ROOT headers
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TSystem.h>


//FastME global manager function
int main(int argc, char *argv[]){
  
  ///Checks for config file
  if(argc < 3 && (argv[1] != help && argv[1] != nc)){
    std::cout<< ansi_red <<"[ERROR]"<< ansi_reset <<" Where's the config file??"<<std::endl;
    return -1;
  }
  
  ///Instantiate the needed variables
  FmeSetup setup; //Class to handle with FastME setups
  TTree *rtree, *ftree;

  ///Takes the config input file and converts it to FastME readable format
  if(argv[1] != help && argv[1] != nc){
    std::cout<<"\n:: ["<<ansi_yellow<<"Your input file: "<< argv[2] <<ansi_reset<<"]"<<std::endl;
    
    //Fills up the 'setup' struct
    ConfigReader((std::string)argv[2], &setup);
    
    ///User can abort analysis if something is wrong in his config file
    std::cout<<ansi_yellow<<" Input file correct?(y/n) "<<ansi_reset;
    std::string aws3;
    std::cin >> aws3;
    if(aws3 == "n") return -1;  
  }

  ///--------------------- Interface manager ------------------------------
  ///To show FastME usage
  if(argv[1] == help)           Helper();
  
  ///To show number of cores available in the local machine
  else if(argv[1] == nc)        FindCores();
  
  ///To format an ntuple in different format
  else if(argv[1] == ff)        FileFormater(setup);
  
  ///To make the Fast Matrix Element analysis and compute discriminant
  else if(argv[1] == fa){
    std::cout<<"\n\n"<<ansi_blue;
    std::cout<<"==============================================================================================="<<std::endl;
    std::cout<<"::::::::::::::::::::::::::[ "<<ansi_cyan<<"Fast Matrix Element Analysis Started"<<ansi_blue<<" ]:::::::::::::::::::::::::::::"<<std::endl;
    std::cout<<"==============================================================================================="<<std::endl;
    std::cout<<ansi_reset;

    ///Calls PhsDrComputer to compute events distance
    rtree = ComputePhsDR(setup);

    ///Calls Discriminator
    ftree = Discriminant(rtree, setup);

    ///Store results (be aware.. the file is handled relative to path where fastme software was called)
    gSystem->Exec("mkdir -p "+setup.OutPath);
    TString resulting_file = setup.OutPath+"/"+setup.OutName+".root";
    TFile *ffile = new TFile(resulting_file,"recreate");
    ftree->Write();
    ffile->Close();
  
    ///-------------------------------------------------------------------------------------------------------------------
    std::cout<<ansi_blue;
    std::cout<<"============================================================================================="<<std::endl;
    std::cout<<":::::::::::::::::::::::::[ "<<ansi_cyan<<"Fast Matrix Element Analysis Finalized"<<ansi_blue" ]::::::::::::::::::::::::::"<<std::endl;
    std::cout<<":: ["<<ansi_yellow<<"Analysis file saved: "<<resulting_file<<ansi_reset<<"]"<<std::endl;
    std::cout<<"============================================================================================="<<std::endl;
    std::cout<<ansi_reset<<"\n\n";
  }

  ///Calls FastME analyzer (plot discriminants, ROC curve, events/discriminant value)
  else if(argv[1] == pr )	StudyResults(setup);
  
  ///Calls event display
  else if(argv[1] == sp)	ShowParticles(setup);
  
  ///Wrong commands gets error and return the helper
  else{
    std::cout<<ansi_red<<"[ERROR]"<<ansi_reset<<" Invalid command '"<<argv[1]<<"'"<<std::endl;
    std::cout<<"These are the available commands:"<<std::endl;
    Helper();
    return -1;
  }
  ///----------------------------------------------------------------------  
 
  //if everything ok, finish well!
  return 0;
}
