///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::[ Fast Matrix Element Executable ]::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



///My Own headers
#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h" //Basic FastME private variables
#include "FastMatrixElement/FastMatrixElement/interface/Cortana.h"     //Config file handler
#include "FastMatrixElement/FastMatrixElement/interface/Librarian.h"      //Handle the original files for reduction and input order indexing
#include "FastMatrixElement/FastMatrixElement/interface/Cartographer.h"   //Mapper real event to MC event
#include "FastMatrixElement/FastMatrixElement/interface/Arbiter.h"        //Discriminator
#include "FastMatrixElement/FastMatrixElement/interface/Composer.h"      //Event generator based on minimum distance method



//ROOT headers
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TSystem.h>


//c++ headers
#include <iostream>
#include <string>



//FastME global manager function
int main(int argc, char *argv[]){

  ///Protection against only software call
  if(argc == 1){
    std::cout<<"\nMissing arguments...\n"<<std::endl;
    Helper();
    return -1;
  }
  ///Protection against wrong command spelling
  if(argv[1] != help && argv[1] != sl && argv[1] != nc && argv[1] != fa && argv[1] != pr && argv[1] != ge){
    std::cout<<"\nBad arguments...\n"<<std::endl;
    Helper();
    return -1;
  }
  
  ///Checks for config file
  if(argc < 3 && (argv[1] != help && argv[1] != nc)){
    std::cout<< ansi_red <<"[ERROR]"<< ansi_reset <<" Where's the config file??"<<std::endl;
    return -1;
  }
  
  ///Instantiate the needed variables
  FmeSetup setup; //Class to handle with FastME setup

  ///Takes the config input file and converts it to FastME readable format
  if(argv[1] != help && argv[1] != nc){
    if(argc < 4){
      //std::cout<<"\n:: ["<<ansi_yellow<<"Your input file: "<< argv[2] <<ansi_reset<<"]"<<std::endl;
    
      //Fills up the 'setup' struct
      ConfigReader((std::string)argv[2], &setup, (std::string)argv[1]);
    
      ///User can abort analysis if something is wrong in his config file
      std::cout<<ansi_yellow<<" Proceed?(y/n) "<<ansi_reset;
      std::string aws3;
      std::cin >> aws3;
      if(aws3 != "y") return -1;  
      
    }
    else ConfigReader((std::string)argv[2], &setup, (std::string)argv[1], (std::string)argv[3]);
  }

  
  ///To show FastME usage
  if(argv[1] == help)           Helper();
  
  
  ///To show number of cores available in the local machine
  else if(argv[1] == nc)        FindCores();
  
  
  ///To make the Fast Matrix Element analysis and compute discriminant
  else if(argv[1] == fa){

    std::cout<<"\n\n"<<ansi_blue;
    std::cout<<"==============================================================================================="<<std::endl;
    std::cout<<"::::::::::::::::::::::::::[ "<<ansi_cyan<<"Fast Matrix Element Analysis Started"<<ansi_blue<<" ]:::::::::::::::::::::::::::::"<<std::endl;
    std::cout<<"==============================================================================================="<<std::endl;
    std::cout<<ansi_reset;

    ///Timming
    TStopwatch t1;
    t1.Start();


    ///Calls PhsDrComputer to compute events distance
    TTree *rtree = Cartographer(setup);


    ///Store results (be aware.. the file is handled relative to path where fastme software is called)
    gSystem->Exec("mkdir -p " + setup.OutPath);
    TString resulting_file = setup.OutPath+"/"+setup.OutName+".root";
    TFile *ffile = new TFile(resulting_file,"recreate");
    rtree->Write();
    ffile->Close();
  
    ///-------------------------------------------------------------------------------------------------------------------
    std::cout<<ansi_blue<<std::endl;
    std::cout<<"==============================================================================================="<<std::endl;
    std::cout<<"::::::::::::::::::::::::::[ "<<ansi_cyan<<"Fast Matrix Element Analysis Finalized"<<ansi_blue<<" ]:::::::::::::::::::::::::::"<<std::endl;
    std::cout<<":: "<<ansi_yellow<<"Analysis file saved: "<<resulting_file<<std::endl;
    std::cout<<ansi_blue<<":: "<<ansi_yellow<< Form("Total Time Analysis: %.2f",t1.RealTime()) <<std::endl;
    std::cout<<ansi_blue<<"==============================================================================================="<<std::endl;
    std::cout<<ansi_reset<<"\n\n";
  }


  ///Calls the Discriminator
  else if (argv[1] == pr){
    Arbiter(setup);
  }
  
  else if (argv[1] == ge){
    Composer(setup);
  }

  
  ///Wrong commands gets error and return the helper
  else{
    std::cout<<ansi_red<<"[ERROR]"<<ansi_reset<<" Invalid command '"<<argv[1]<<"'"<<std::endl;
    std::cout<<"These are the available commands:"<<std::endl;
    Helper();
    return -1;
  }

 
  //if everything ok, finish well!
  return 0;
}
