///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::: FastME Definitions for Interface :::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#include <TString.h>
#include <string>
#include <vector>

///-------------  Define the key-words  ------------------
vector<string> fme_keywords = {
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
  "verbose_level"
};

vector<int> ksize = {
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
  14
};

TString Data_Path;

vector<string> MCs;
///-------------------------------------------------------
