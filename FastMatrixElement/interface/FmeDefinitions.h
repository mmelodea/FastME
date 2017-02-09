///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::[ FastME Global Definitions ]:::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef FmeDefinitions_h
#define FmeDefinitions_h

#include <TFile.h>
#include <TString.h>
#include <string>
#include <vector>



///System output colors
static std::string ansi_red    = "\x1b[31m";
static std::string ansi_green  = "\x1b[32m";
static std::string ansi_yellow = "\x1b[33m";
static std::string ansi_blue   = "\x1b[34m";
static std::string ansi_violet = "\x1b[35m";
static std::string ansi_cyan   = "\x1b[36m";
static std::string ansi_reset  = "\x1b[0m";




///------- Define global struct to store user setup ----------
struct FmeSetup{
  std::vector<std::string>	vDatas;
  Int_t				NData;
  TString			TTreeName;
  TString			IdBranch;
  TString			PtBranch;
  TString			EtaBranch;
  TString			PhiBranch;
  std::vector<std::string>	vMCs;
  Int_t				NCores;
  TString			PhSDrMethod;
  TString			SetFlavorConstraint;
  TString			OutName;
  TString			OutPath;
  Int_t				DTLimit;
  Float_t			MCLimit;
  Double_t			ScaledPt;
  Double_t			ScaledEta;
  Double_t			ScaledPhi;
  TString			FmeFile;
  std::vector<int>		SigData;
  std::vector<int>              BkgData;
  std::vector<int>              SigMC;
  std::vector<int>              BkgMC;
  TString			XSScale;
  std::vector<float>		DataXS;
  std::vector<float>		McXS;
  Int_t				GenNEv;
  Double_t			SDrCondition;
  Double_t                      BDrCondition;
  Int_t 			MaxGenTrials;
  TString			OutGenName;
  Int_t				Verbose = 1;
};



///-------------  Define global key-words  ------------------
static std::vector<std::string> fme_keywords = {
  "data_path",
  "mc_path",
  "tree_name",
  "id_branch_name",
  "pt_branch_name",
  "eta_branch_name",
  "phi_branch_name",
  "outfile_name",
  "outfile_path",
  "phs_dr_method",
  "flavor_constraint",
  "n_cores",
  "data_limit",
  "mc_limit",
  "scale_dPt",
  "scale_dEta",
  "scale_dPhi",
  "fme_file",
  "sig_data",
  "bkg_data",
  "sig_mc",
  "bkg_mc",
  "xs_scale",
  "data_xs",
  "mc_xs",
  "gen_nevents",
  "sdr_condition",
  "bdr_condition",
  "max_gen_trials",
  "out_gen_name",
  "verbose_level"
};



static std::vector<int> ksize = {
  10,	///data_path
  8,	///mc_path
  10,	///tree_name
  15,	///id_branch_name
  15,	///pt_branch_name
  16,	///eta_branch_name
  16,	///phi_branch_name
  13,	///outfile_name
  13,	///outfile_path
  14,	///phs_dr_method
  18,	///flavor_constraint
  8,	///n_cores
  11,	///data_limit
  9,	///mc_limit
  10,	///scale_dPt
  11,	///scale_dEta
  11,	///scale_dPhi
  9,	///fme_file
  9,	///sig_data
  9,	///bkg_data
  7,	///sig_mc
  7,	///bkg_mc
  9,	///xs_scale
  8,	///data_xs
  6,	///mc_xs
  12,	///gen_nevents
  14,	///sdr_condition
  14,   ///bdr_condition
  15,	///max_gen_trials
  13,	///out_gen_name
  14	///verbose_level
};


#endif
