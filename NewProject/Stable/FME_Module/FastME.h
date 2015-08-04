#ifndef FME_h
#define FME_h

#include <TTree.h>
#include <TString.h>

class FME{
  
  public:
    FME(){}
    ~FME(){}
    int launchFME(TString Final_State, TString Model, TString Out_Name, TString Data_Path, TString MC_Sig_Path, TString MC_Bkg_Path, TString Tree_Name, TString Branch_Name, TString Resolution);
  
  private:
    int FS_4l(TString Model, TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name, TString Resolution);
    int FS_4l2j(TString Model, TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name, TString Resolution);
    int FS_lv2j(TString Model, TString Out_Name, TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name, TString Resolution);    

};

#endif