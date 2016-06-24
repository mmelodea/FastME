///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Converter to FastME Input File Format ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef FileFormater_h
#define FileFormater_h

#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <string>

///Tool to convert any root file to the FastME format
void FileFormater(FmeSetup UserSetup){
  //put option to use tproc.. so many files can be converted at the same time
  std::cout<<"Option under construction yet..."<<std::endl;  
  
  //Get how many files to be converted
  //Int_t nCnvFiles = UserSetup.CnvFile.size();
  
  //for(Int_t i_f=0; i_f<nCnvFiles; i_f++){
    //TFile *infile = TFile::Open(UserSetup.CnvFile);
    //TTree *intree = (TTree*)infile->Get(CnvTree);
    //intree->SetBranchAddress(UserSetup.CnvBranch[],&);
  //}
  
  return;
}

#endif