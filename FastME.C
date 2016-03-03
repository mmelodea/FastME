///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::: FastME Interface :::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#include "FmeInterfaceDef.h"
#include "ComputePhsDR.C"
#include "Discriminant.C"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <string>

int FastME(string config){
  
  ///_____________________________ For timming the process ________________________________________________________
  cout<<"\n\n";
  cout<<"=============================================================================================="<<endl;
  cout<<"::::::::::::::::::::::::::[ Fast Matrix Element Analysis Started ]::::::::::::::::::::::::::::"<<endl;
  cout<<"=============================================================================================="<<endl;
  cout<<" Certify RAM memory available? (y/n) ";
  string aws1;
  cin >> aws1;
  if(aws1 == "y"){
    cout<<" [RAM memory info]"<<endl;
    gSystem->Exec("free -m");
    cout<<" Proceed? (y/n) ";
    string aws2;
    cin >> aws2; 
    if(aws2 != "y"){
      cout<<" Stoped due to insuficient RAM memory availability!"<<endl;
      return -1;
    }
  }
  cout<<"\n:: Your input file..."<<endl;
  ///______________________________________________________________________________________________________________
  
  ///Define variables used in the analysis
  MCs.clear();
  TString PhSDr_Method, FlavorConstraint, out_name, out_path;
  TString TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch;
  vector<TString> MC_Names;
  Int_t N_FSParticles= 4, verbose= 1, DT_Limit= -1;
  Float_t MC_Limit= -1, scale_dPt = 50., scale_dEta = 5.;
  UInt_t N_Cores= 1;

  ///______________ Extract the need info from txt file _______________________________________________________________
  ///Opens the config file to get the user configuration
  ifstream inFile(config.c_str());
  if(!inFile){
    cout<<"Error in file!";
    return -1;
  }
  int nkeys=0;
  string line;
  while( getline(inFile,line) ){
    if(line.find("#") != string::npos) continue;
    for(int k=0; k<(int)fme_keywords.size(); k++){
      if(line.find(fme_keywords[k]) != string::npos){
	nkeys++;
        line.erase(line.begin(),line.begin()+ksize[k]);
        if(k==0 || k==1 || k==2 || k==3 || k==8 || k==9 || k==13 || k==14 || k==15 || k==16 || k==17)
	  cout <<":: "<< fme_keywords[k] <<"\t\t\t\t"<< line << endl;
	else cout <<":: "<< fme_keywords[k] <<"\t\t\t"<< line << endl;
        if(fme_keywords[k] == "data_path") Data_Path = line;
	if(fme_keywords[k] == "mc_path") MCs.push_back(line);
	if(fme_keywords[k] == "mc_name") MC_Names.push_back(line);
	if(fme_keywords[k] == "tree_name") TreeName = line;
	if(fme_keywords[k] == "mc_type_branch_name") McType_branch = line;
	if(fme_keywords[k] == "id_branch_name") Id_branch = line;
	if(fme_keywords[k] == "pt_branch_name") Pt_branch = line;
	if(fme_keywords[k] == "eta_branch_name") Eta_branch = line;
	if(fme_keywords[k] == "outfile_name") out_name = line;
	if(fme_keywords[k] == "outfile_path") out_path = line;
	if(fme_keywords[k] == "phs_dr_method") PhSDr_Method = line;
	if(fme_keywords[k] == "n_fs_particles") N_FSParticles = stoi(line);
	if(fme_keywords[k] == "flavor_constraint") FlavorConstraint = line;
	if(fme_keywords[k] == "n_cores") N_Cores = stoi(line);
	if(fme_keywords[k] == "data_limit") DT_Limit = stoi(line);
	if(fme_keywords[k] == "mc_limit") MC_Limit = stof(line);
        if(fme_keywords[k] == "scale_dPt") scale_dPt = stof(line);
        if(fme_keywords[k] == "scale_dEta") scale_dEta = stof(line);
	if(fme_keywords[k] == "verbose_level") verbose = stoi(line);
      }
    }
  }
  if(nkeys < ((int)fme_keywords.size()-2)){
    cout<<"Missing key-word! Check your input file!";
    return -1;
  }
  cout<<" Input file correct?(y/n) ";
  string aws3;
  cin >> aws3;
  if(aws3 == "n") return -1;
  ///__________________________________________________________________________________________________________________
  
  ///Getting some numbers
  TFile *fData = TFile::Open(Data_Path);
  TTreeReader tmpReader(TreeName,fData);
  Int_t nData = tmpReader.GetEntries(true);
  if(DT_Limit != -1 && DT_Limit <= nData)
    nData = DT_Limit;

  const Int_t N_MCT = MC_Names.size();
  const Int_t N_MC = MCs.size();
  Int_t NMCEV[N_MC];
  for(Int_t ne=0; ne<N_MC; ne++){
    TFile *tmpf = TFile::Open((TString)MCs[ne]);
    TTreeReader tmpReader(TreeName,tmpf);
    if(MC_Limit == -1)
      NMCEV[ne] = tmpReader.GetEntries(true);
    if(MC_Limit != -1 && MC_Limit < 1)
      NMCEV[ne] = (Int_t)(MC_Limit*tmpReader.GetEntries(true));
    if(MC_Limit != -1 && MC_Limit >= 1)
      NMCEV[ne] = MC_Limit;
    tmpf->Close();
  }

  cout<<"______________________________________________________________________________________________"<<endl;
  cout<< Form(":: Data Events:      %i",nData) <<endl;
  cout<< Form(":: MC Samples:       %i\t[",N_MC); 
            for(Int_t ne=0; ne<N_MC; ne++){
              if( ne<(N_MC-1) ) cout<< ne << "= " << NMCEV[ne] << ",  ";
              if( ne==(N_MC-1) ) cout<< ne << "= " << NMCEV[ne] << "]" <<endl;
            }
  cout<< Form(":: Final State:      %i  Objects",N_FSParticles) <<endl;
  cout<< Form(":: Cores to Use:     %i  Cores",N_Cores) <<endl;
  cout<<"----------------------------------------------------------------------------------------------"<<endl;

  if(aws1 == "y"){
    cout<<" [RAM info after full compilation]:"<<endl;
    gSystem->Exec("free -m");
    cout<<" Continue? (y/n) ";
    string aws4;
    cin >> aws4;
    if(aws4 != "y"){
      cout<<" Stoped stopped by user!"<<endl;
      return -1;
    }
  }
  cout<<" [Analysing events...]"<<endl;
  ///--------------------------------------------------------------------------------------------------------------
  
  ///Calls TProcPool to compute distances between events
  TTree *rtree = 
  ComputePhsDR(fData, nData, TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch, MCs, N_MCT, N_Cores,
	       N_FSParticles, PhSDr_Method, FlavorConstraint, MC_Limit, scale_dPt, scale_dEta, verbose);
  
  ///Calls discriminator
  TTree *ftree = Discriminant(rtree, N_Cores, nData, N_MCT, verbose);
  
  ///Saving FastME results
  gSystem->Exec("mkdir -p "+out_path);
  TFile *tmp = TFile::Open(out_path+"/"+out_name+".root","recreate");
  ftree->Write();
  tmp->Close();
  
  cout<<"=============================================================================================="<<endl;
  cout<<":::::::::::::::::::::::::[ Fast Matrix Element Analysis Finalized ]:::::::::::::::::::::::::::"<<endl;
  cout<<"=============================================================================================="<<endl;

  return 0;
}
