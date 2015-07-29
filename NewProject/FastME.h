#ifndef FME_h
#define FME_h

#include <TTree.h>
#include <TString.h>

class FME{
  
  public:
    FME(){}
    ~FME(){}
    TTree* launch_FME(TString type, TString out_name, TTree *Data_tree, TTree *MC_Sig_tree, TTree *MC_Bkg_tree, TString Branch_Name);
  
  private:
    TTree* FS4l_DR(TString out_name, TTree *Data_tree, TTree *MC_Sig_tree, TTree *MC_Bkg_tree, TString Branch_Name);
    //TTree* FS4l2j_DR(string out_name, TTree *Data_tree, TTree *MC_Sig_tree, TTree *MC_Bkg_tree, string Branch_Name);
    //TTree* FSlv2j_DR(string out_name, TTree *Data_tree, TTree *MC_Sig_tree, TTree *MC_Bkg_tree, string Branch_Name);

};

#endif