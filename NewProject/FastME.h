#ifndef FME_h
#define FME_h

#include <TTree.h>
#include <TString.h>

//private:
  //double std_DR(const int comp, Float_t Data[comp][3][2], Float_t MC[comp][3][2]);

class FME{
  
  public:
    FME(){}
    ~FME(){}
    int launchFME(TString Final_State, TString Out_Name, TString Data_Path, TString MC_Sig_Path, TString MC_Bkg_Path, TString Tree_Name, TString Branch_Name);
  
  private:
    int FS4l(TString Out_Name,TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
    int FS4l2j(TString Out_Name,TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
    int FSlv2j(TString Out_Name,TTree *Data_Tree, TTree *MC_Sig_Tree, TTree *MC_Bkg_Tree, TString Branch_Name);
    
};

#endif