///####################################################################################################################
///=========================================[ Fast Matrix Element Module ]=============================================
///=====================================[ Code Author: Miqueias M. de Almeida ]========================================
///####################################################################################################################



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TObjArray.h>

///Headers to TProcPool
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TProcPool.h>
#include <PoolUtils.h>
#include <TSystem.h>


#define pi		3.14159265358979312
///Scale Factors to normalize deltas
#define scale_dPt	50.
#define scale_dEta	5.
#define scale_dPhi	pi

using namespace std;

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
  14
};
vector<string> null = {""};
///-------------------------------------------------------



///========== Compute Discriminant Value ===============
///Based on Distance
Double_t PsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}
///=====================================================


///######################################## Main FastME Function ######################################################
int FastME(TString Data_Path="", vector<string> MCs=null){

  
  ///_____________________________ For timming the process ________________________________________________________
  TStopwatch t;
  t.Start();
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
  if(MCs == null) MCs.clear();
  TString PhSDr_Method, FlavorConstraint, out_name, out_path;
  TString TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch;
  vector<TString> MC_Names;
  Int_t N_FSParticles= 4, verbose= 1, DT_Limit= -1;
  Float_t MC_Limit= -1;
  UInt_t N_Cores= 1;

  ///______________ Extract the need info from txt file _______________________________________________________________
  ///Opens the config file to get the user configuration
  string config = "fme_config.dat";
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
	cout <<":: "<< fme_keywords[k] <<"\t\t"<< line << endl;
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
  
  
  ///TProcPool declaration to objects to be analised  
  auto workItem = [fData, nData, TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch, N_MCT, N_FSParticles,
		   PhSDr_Method, FlavorConstraint, MC_Limit, verbose](TTreeReader &tread) -> TObject* {
    TStopwatch t2;
        
    ///Addresses the MC branches to be used
    TTreeReaderValue<Int_t>    McType(tread, McType_branch); ///McType for Signal=0 and Background >0
    TTreeReaderArray<Int_t>    McId(tread, Id_branch);
    TTreeReaderArray<Double_t> McPt(tread, Pt_branch);
    TTreeReaderArray<Double_t> McEta(tread, Eta_branch);
    
    ///Addresses the Data branches to be used
    TTreeReader refReader(TreeName,fData);
    TTreeReaderArray<Int_t>    DataId(refReader, Id_branch);
    TTreeReaderArray<Double_t> DataPt(refReader, Pt_branch);
    TTreeReaderArray<Double_t> DataEta(refReader, Eta_branch);


    ///Tree to store the results from analysis
    Int_t iEvent, TMcType, Indice;
    Double_t Mdist;
    TTree *fme_tree = new TTree("fme_tree","temporary");
    fme_tree->SetDirectory(0);
    fme_tree->Branch("iEvent",&iEvent,"iEvent/I");
    fme_tree->Branch("Mdist",&Mdist,"Mdist/D");
    fme_tree->Branch("TMcType",&TMcType,"TMcType/I");
    fme_tree->Branch("Indice",&Indice,"Indice/I");
    
    ///Loop on Data events
    for(Int_t dt=0; dt<nData; dt++){
      if( verbose != 0 && ((dt!= 0 && dt%(nData/10) == 0) || (nData-dt) == 1) ){ 
	cout<< Form(":: [Remaining]:  %i Events\t\t[Elapsed]:  ",nData-dt);
	t2.Stop();
	t2.Print();
	t2.Continue();
      }
      refReader.SetEntry(dt); ///Move on Data loop
      Double_t min_distance_Min = 1.e15;
      Double_t min_distance_Med = 1.e15;
      Int_t imc_min = -1;
      Int_t f_type=-99;
      Int_t nMonteCarlo = tread.GetEntries(true);
      if(MC_Limit != -1 && MC_Limit >= 1) nMonteCarlo = MC_Limit;
      if(MC_Limit != -1 && MC_Limit < 1) nMonteCarlo = (Int_t)(MC_Limit*nMonteCarlo);
      
      for(Int_t mc=0; mc<nMonteCarlo; mc++){
	tread.SetEntry(mc); ///Move on MC loop
        bool acept = true;
	
      ///==============================================================================================================
      ///::::::::::::::::::::::::: Fast Matrix Element methods to compute Data - MC distance ::::::::::::::::::::::::::
      ///::::::::::::::::::::::::::::::::::::: Currently 2 Methods Available ::::::::::::::::::::::::::::::::::::::::::
      ///==============================================================================================================
	Double_t event_distance_Min= -99, event_distance_Med= -99, event_distance= -99;
	Double_t SumMed_dPt2 = 0, SumMed_dEta2 = 0;
	Double_t SumMin_dPt2 = 0, SumMin_dEta2 = 0;
	
	//stores flags to sinalize when a MC object is already selected
	vector<int> vmin_imc;
	for(int sl=0; sl<N_FSParticles; sl++) vmin_imc.push_back(-1);
	
	for(int idt=0; idt<N_FSParticles; idt++){
	  Double_t min_particles_distance = 1.E15;
	  Double_t particles_distance = -1.;
	  int min_imc = -1;
    
	  Int_t nsame_flavor = 0;
	  Double_t tmp_dPt = 0, tmp_dEta = 0;
	  for(int imc=0; imc<N_FSParticles; imc++){
	    ///Avoid different Data-MC particles comparison
	    if(FlavorConstraint == "true" && DataId[idt] != McId[imc]) continue;
	    ///Avoid leptons-jets comparison
	    else if(FlavorConstraint == "false"){
	      if(abs(DataId[idt])== 11 && (abs(McId[imc])!= 11 || abs(McId[imc])!= 13)) continue;
	      if(abs(DataId[idt])== 13 && (abs(McId[imc])!= 11 || abs(McId[imc])!= 13)) continue;
	      if(abs(McId[idt])== 11 && (abs(DataId[imc])!= 11 || abs(DataId[imc])!= 13)) continue;
	      if(abs(McId[idt])== 13 && (abs(DataId[imc])!= 11 || abs(DataId[imc])!= 13)) continue;
	    }
	    ///Compute preliminary particles distance
	    Double_t dPt  = (DataPt[idt]-McPt[imc])/(scale_dPt);
	    Double_t dEta = (DataEta[idt]-McEta[imc])/(scale_dEta);
	

	  ///_______________________ For proximity comparison method __________________________________________________
	    if( PhSDr_Method == "mindr")
	      if(vmin_imc[imc] == -1){
		particles_distance = sqrt(dPt*dPt + dEta*dEta);
		if( verbose == 3 ) cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<"   part_dist: "<<particles_distance<<endl;
		if(particles_distance < min_particles_distance){
		  min_imc = imc;
		  min_particles_distance = particles_distance;
		}
	      }
	  ///¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
	  
	  
  	  ///________________________ Only for media comparison method ________________________________________________
	    if(PhSDr_Method == "media"){
	      if(nsame_flavor == 0){
		tmp_dPt  = dPt/scale_dPt;
		tmp_dEta = dEta/scale_dEta;
	      }
	      else{
		///Repair the previous one
		SumMed_dPt2  += pow(0.5*tmp_dPt,2);
		SumMed_dEta2 += pow(0.5*tmp_dEta,2);
		///Append the new one
		SumMed_dPt2  += pow(0.5*tmp_dPt,2);
		SumMed_dEta2 += pow(0.5*tmp_dEta,2);
	      }
	      nsame_flavor++;
	      if( verbose == 3 )
		cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<endl;
	    }
	  ///¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨

	  }///Ends MC event loop
	  
	  if(PhSDr_Method == "mindr"){
	    ///Monitor of chosen MCs to avoid object recounting and wrong pairing
	    if(min_imc == -1){
	      acept = false;
	      break;
	    }
	    vmin_imc[min_imc] = 1;//changes the flag for current MC object
	    if( verbose == 3 ) cout<<"Chosen->>  MCPos: "<<min_imc<<"   ID: "<<McId[min_imc]<<endl;
	    ///For proximity comparison method
	    SumMin_dPt2 += pow( (DataPt[idt]-McPt[min_imc])/(scale_dPt), 2 );
	    SumMin_dEta2 += pow( (DataEta[idt]-McEta[min_imc])/(scale_dEta), 2 );
	  }
	  
	  if(PhSDr_Method == "media" && nsame_flavor == 1){
	    SumMed_dPt2  += tmp_dPt*tmp_dPt;
	    SumMed_dEta2 += tmp_dEta*tmp_dEta;
	  }
	}///Ends Data event loop
	
	///Compute final Data-MC events distance & searches for minimum distance
	if(PhSDr_Method == "mindr" && acept == true){
	  event_distance_Min = sqrt(SumMin_dPt2 + SumMin_dEta2);
	  if(event_distance_Min < min_distance_Min){
	    min_distance_Min = event_distance_Min;
	    imc_min = mc;
	  }
	  if( verbose > 2 ) cout<<"Event_distance(MinDr) = "<<event_distance_Min<<endl;  
	}
	
	if(PhSDr_Method == "media" && SumMed_dPt2 > 0){
	  event_distance_Med = sqrt(SumMed_dPt2 + SumMed_dEta2);
	  if(event_distance_Med < min_distance_Med){
	    min_distance_Med = event_distance_Med;
	    imc_min = mc;
	  }
	  if( verbose > 2 ) cout<<"Event_distance (Media) = "<<event_distance_Med<<endl;  
	}
	///Stores the MC type
	f_type = *McType;
      ///================================================================================================================
      ///================================================================================================================
	
      }///End MC sample loop

      
      ///Stores the minimum distances found
      iEvent = dt;
      if(PhSDr_Method == "mindr"){
	Mdist = min_distance_Min;
	TMcType = f_type;
	Indice  = imc_min;
	if( verbose > 1 ) cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Min<<endl;
      }
      if(PhSDr_Method == "media"){
	Mdist = min_distance_Med;
	TMcType = f_type;
	Indice  = imc_min;
	if( verbose > 1 ) cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance("<<imc_min<<"): "<<min_distance_Med<<endl;
      }
      
      fme_tree->Fill();
    }///End Data sample loop
    
    t2.Stop();
    delete fData;
    return fme_tree;
  };
  
  ///Calls analysis through TProcPool
  TProcPool workers(N_Cores);
  TTree *mtree = (TTree*)workers.ProcTree(MCs, workItem);
  Int_t iEvent, TMcType, Indice;
  Double_t Mdist;
  mtree->SetBranchAddress("iEvent",&iEvent);
  mtree->SetBranchAddress("Mdist",&Mdist);
  mtree->SetBranchAddress("TMcType",&TMcType);
  mtree->SetBranchAddress("Indice",&Indice);
  Int_t fentries = mtree->GetEntries();

  ///_______________________ Compute discriminant from MDMCED _____________________________________________________
  cout<<":: [Distance Computing Time]: "; t.Stop(); t.Print(); t.Continue();
  cout<<"\n::::::::::::::::::::::::::::::::[ Computing discriminant ]::::::::::::::::::::::::::::::::::::"<<endl;
  ///--------------------------------------------------------------------------------------------------------------

  vector<Int_t> McIndex, McCat;
  Double_t Global_PsbDist;
  vector<Double_t> MinDist, Local_PsbDist;
  TTree *ftree = new TTree("FastME","Fast Matrix Element Analysis Results");
  ftree->SetDirectory(0);
  ftree->Branch("McIndex","vector<Int_t>",&McIndex);
  ftree->Branch("McCat","vector<Int_t>",&McCat);
  ftree->Branch("MinDist","vector<Double_t>",&MinDist);
  ftree->Branch("Global_PsbDist",&Global_PsbDist,"Global_PsbDist/D");
  ftree->Branch("Local_PsbDist","vector<Double_t> Local_PsbDist",&Local_PsbDist);

  ///Find the tree sectors
  Int_t TreeSectors[N_Cores];
  if(fentries % N_Cores != 0){
    cout<<"[Error] Something gone wrong... non-integer tree sectors!!"<<endl;
    throw exception();
  }
  
  Int_t EndSector = fentries/N_Cores; //How many entries in each core job                                     
  for(Int_t ic=0; ic<(Int_t)N_Cores; ic++) TreeSectors[ic] = ic*EndSector;
  
  ///Getting results from analysis
  for(Int_t data=0; data<nData; data++){
    if( verbose != 0 && data%(nData/10) == 0)
      cout<< Form(":: [Remaining]:   %i Events",nData-data) <<endl;

    Int_t MinSigIndex = -99, GMinBkgIndex = -99;
    Double_t min_dr_sig = 1.e15, global_min_dr_bkg = 1.e15;
    Double_t local_min_dr_bkg[N_MCT], local_min_bkg_index[N_MCT];
    McIndex.clear();
    McCat.clear();
    MinDist.clear();
    Local_PsbDist.clear();
    for(Int_t p=0; p<N_MCT; p++){
      local_min_bkg_index[p] = -99;
      local_min_dr_bkg[p] = 1.e15;
      McCat.push_back( -99 );
      McIndex.push_back( -99 );
      MinDist.push_back( -99. );
      Local_PsbDist.push_back( -99. );
    }
    for(Int_t ic=0; ic<(Int_t)N_Cores; ic++){
      mtree->GetEntry(TreeSectors[ic]+data);//TreeSectors aligns the results from different cores
      if(iEvent != data){
	cout<<"[Error] Something gone wrong... iEvent != data"<<endl;
	throw exception();
      }
      
      ///Finds closet MC Signal
      if(TMcType == 0)
	if( Mdist < min_dr_sig ){
	  min_dr_sig = Mdist;
	  MinSigIndex = Indice;
	  McCat[0] = 0;
	  //cout<<"Indice= "<<Indice<<endl;
	}

      ///Finds closet MC Background
      if(TMcType > 0){
	///The general most close MC Background
        if( Mdist < global_min_dr_bkg ) global_min_dr_bkg = Mdist;
	///Each MC Background
	if( Mdist < local_min_dr_bkg[TMcType-1] ){
	  McCat[TMcType] = TMcType;
	  local_min_dr_bkg[TMcType-1] = Mdist;
	  local_min_bkg_index[TMcType-1] = Indice;
	}
      }
    }

    MinDist[0] = min_dr_sig;
    McIndex[0] = MinSigIndex;
    Global_PsbDist = PsbD(min_dr_sig, global_min_dr_bkg);

    for(Int_t im=0; im<N_MCT; im++){
      McIndex[im+1] = local_min_bkg_index[im];
      MinDist[im+1] = local_min_dr_bkg[im];
      Local_PsbDist[im+1] = PsbD(min_dr_sig, local_min_dr_bkg[im]);
    }
    if( verbose > 1 )
      cout<< Form("GSigMin:   %f\t\tGBkgMin:   %f\t\tGPsbDMinDist:   %f", min_dr_sig, global_min_dr_bkg, Global_PsbDist) << endl;
    
    
    ftree->Fill();
  }

  ///________________________________ Stoping timming ________________________________________________________
  cout<<"\n::::::::::::::::::::::::::::::::::::[ Process Finished ]::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: [Analysis Total Time]: "; t.Stop(); t.Print();
  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"<<endl;
  ///---------------------------------------------------------------------------------------------------------

  ///Saving FastME results
  gSystem->Exec("mkdir -p "+out_path);
  TFile *tmp = TFile::Open(out_path+"/"+out_name+".root","recreate");
  ftree->Write();
  tmp->Close();
  
  return 0;
}
///********************************************************************************************************************
