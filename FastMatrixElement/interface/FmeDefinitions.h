///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::: FastME Definitions for Interface :::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifndef FmeDefinitions_h
#define FmeDefinitions_h

#include <TFile.h>
#include <TString.h>
#include <string>
#include <vector>



///------- Define global struct to store user setup ----------
struct FmeSetup{
  TFile				*DataFile;
  Int_t				NData;
  TString			TTreeName;
  TString			McTypeBranch;
  TString			IdBranch;
  TString			PtBranch;
  TString			EtaBranch;
  std::vector<std::string>	vMCs;
  Int_t				NMCT;
  UInt_t			NCores;
  Int_t				NFSParticles;
  TString			PhSDrMethod;
  TString			SetFlavorConstraint;
  TString			OutName;
  TString			OutPath;
  Int_t				DTLimit;
  Float_t			MCLimit;
  Double_t			Scale_dPt;
  Double_t			Scale_dEta;
  TString			ScaleMethod;
  Int_t				Verbose = 1;
};


///-------------  Define global key-words  ------------------
static std::vector<std::string> fme_keywords = {
  "data_path",
  "mc_path",
  "mc_name",
  "tree_name",
  "mc_type_branch_name",
  "id_branch_name",
  "pt_branch_name",
  "eta_branch_name",
  "outfile_name",
  "outfile_path",
  "phs_dr_method",
  "n_fs_particles",
  "flavor_constraint",
  "n_cores",
  "data_limit",
  "mc_limit",
  "scale_dPt",
  "scale_dEta",
  "scale_method"
  "verbose_level"
};

static std::vector<int> ksize = {
  10,
  8,
  8,
  10,
  20,
  15,
  15,
  16,
  13,
  13,
  14,
  15,
  18,
  8,
  11,
  9,
  10,
  11,
  13,
  14
};


#endif
