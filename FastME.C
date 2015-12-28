///####################################################################################################################
///=========================================[ Fast Matrix Element Module ]=============================================
///=====================================[ Code Author: Miqueias M. de Almeida ]========================================
///####################################################################################################################



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStopwatch.h>

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
  cout<<" [RAM memory info]"<<endl;
  gSystem->Exec("free -m");
  cout<<" Certify RAM memory available? (y/n)";
  string aws;
  cin >> aws; 
  if(aws != "y"){
    cout<<" Stoped due to insuficient RAM memory availability!"<<endl;
    return -1;
  }
  cout<<"\n:: Your input file..."<<endl;
  ///______________________________________________________________________________________________________________
  
  ///Define variables used in the analysis
  if(MCs == null) MCs.clear();
  TString PhSDr_Method, FlavorConstraint, out_name, out_path;
  TString TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch;
  vector<TString> MC_Names;
  Int_t N_FSParticles= 4, verbose= 1;
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
        if(fme_keywords[k] == "mc_limit") MC_Limit = stof(line);
	if(fme_keywords[k] == "verbose_level") verbose = stoi(line);
      }
    }
  }
  if(nkeys < ((int)fme_keywords.size()-2)){
    cout<<"Missing key-word! Check your input file!";
    return -1;
  }
  ///__________________________________________________________________________________________________________________
  
  ///Getting some numbers
  TFile *fData = TFile::Open(Data_Path);
  TTreeReader tmpReader(TreeName,fData);
  Int_t nData = tmpReader.GetEntries(true);

  Int_t N_MCT = MC_Names.size();
  Int_t N_MC = MCs.size();
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
  cout<<" [RAM info]:"<<endl;
  gSystem->Exec("free -m");
  cout<<" Continue proccess? (y/n)";
  string aws2;
  cin >> aws2;
  if(aws2 != "y"){
    cout<<" Stoped stopped by user!"<<endl;
    return -1;
  }
  cout<<" [Analysing events...]"<<endl;  
  ///--------------------------------------------------------------------------------------------------------------
  
  
  ///TProcPool declaration to objects to be analised  
  auto workItem = [fData, nData, TreeName, McType_branch, Id_branch, Pt_branch, Eta_branch, N_MCT, N_FSParticles,
		   PhSDr_Method, FlavorConstraint, MC_Limit, verbose](TTreeReader &tread) -> TObject* {
    TStopwatch t2;
        
    ///Defines 2D histograms to stores minimum distances
    TH2D *mdists = new TH2D("mdists","Minimum Data-MC distances found",nData,0,nData,N_MCT,0,N_MCT);
    mdists->SetDirectory(0);    

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


    ///Loop on Data events
    for(Int_t dt=0; dt<nData; dt++){
      if( verbose != 0 && ((dt!= 0 && dt%(nData/10) == 0) || (nData-dt) == 1) ){ 
	cout<< Form(":: [Remaining]:  %i Events\t\t[Elapsed]:  ",nData-dt);
	t2.Stop();
	t2.Print();
	t2.Continue();
      }
      refReader.SetEntry(dt); ///Move on Data loop
      Double_t min_distance_Min = 1.E15;
      Double_t min_distance_Med = 1.E15;
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
	vector<Double_t> SumMed_dPt, SumMed_dEta, SumMin_dPt, SumMin_dEta;
	vector<int> vmin_imc;
	SumMin_dPt.clear();
	SumMin_dEta.clear();
	SumMed_dPt.clear();
	SumMed_dEta.clear();
	vmin_imc.clear();
	
	for(int idt=0; idt<N_FSParticles; idt++){
	  Double_t min_particles_distance = 1.E15;
	  Double_t particles_distance = -1.;
	  int min_imc = -1;
    
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
	    Double_t dPt  = (DataPt[idt]-McPt[imc])/(scale_dPt);
	    Double_t dEta = (DataEta[idt]-McEta[imc])/(scale_dEta);
	

	  ///_______________________ For proximity comparison method __________________________________________________
	    if( PhSDr_Method == "mindr"){
	      bool mc_approved = true;
	      if(int(vmin_imc.size()) > 0)
		for(int g=0; g<int(vmin_imc.size()); g++)
		  if(imc == vmin_imc[g]) mc_approved = false;
	
	      if(mc_approved == true){
		particles_distance = sqrt(dPt*dPt + dEta*dEta);
		if( verbose == 3 ) cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<"   part_dist: "<<particles_distance<<endl;
		if(particles_distance < min_particles_distance){
		  min_imc = imc;
		  min_particles_distance = particles_distance;
		}
	      }
	    }
	  ///¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
	  
	  
  	  ///________________________ Only for media comparison method ________________________________________________
	    if(PhSDr_Method == "media"){
	      SumMed_dPt.push_back(dPt/2.);
	      SumMed_dEta.push_back(dEta/2.);
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
	    vmin_imc.push_back(min_imc);
	    if( verbose == 3 ) cout<<"Chosen->>  MCPos: "<<min_imc<<"   ID: "<<McId[min_imc]<<endl;
	    ///For proximity comparison method
	    SumMin_dPt.push_back( (DataPt[idt]-McPt[min_imc])/(scale_dPt) );
	    SumMin_dEta.push_back( (DataEta[idt]-McEta[min_imc])/(scale_dEta) );
	  }
	}///Ends Data event loop  
	
	///Compute final Data-MC events distance & searches for minimum distance
	Double_t sum_dPt2 = 0, sum_dEta2 = 0;	
	if(PhSDr_Method == "mindr" && acept == true){
	  for(int n=0; n<(int)SumMin_dPt.size(); n++){
	    sum_dPt2  += SumMin_dPt[n]*SumMin_dPt[n];
	    sum_dEta2 += SumMin_dEta[n]*SumMin_dEta[n];
	  }	  
	  event_distance_Min = sqrt(sum_dPt2 + sum_dEta2);
	  if(event_distance_Min < min_distance_Min) min_distance_Min = event_distance_Min;
	  if( verbose > 1 ) cout<<"Event_distance(MinDr) = "<<event_distance_Min<<endl;  
	}
	
	if(PhSDr_Method == "media" && SumMed_dPt.size() > 0){
	  sum_dPt2 = 0; sum_dEta2 = 0; 
	  for(int n=0; n<(int)SumMed_dPt.size(); n++){
	    sum_dPt2  += SumMed_dPt[n]*SumMed_dPt[n];
	    sum_dEta2 += SumMed_dEta[n]*SumMed_dEta[n];
	  }	  
	  event_distance_Med = sqrt(sum_dPt2 + sum_dEta2);
	  if(event_distance_Med < min_distance_Med) min_distance_Med = event_distance_Med;
	  if( verbose > 1 ) cout<<"Event_distance (Media) = "<<event_distance_Med<<endl;  
	}
	///Stores the MC type
	f_type = *McType;
      ///================================================================================================================
      ///================================================================================================================
	
      }///End MC sample loop

      
      ///Stores the minimum distances found
      if(PhSDr_Method == "mindr") mdists->Fill(dt,f_type,min_distance_Min);
      if(PhSDr_Method == "media") mdists->Fill(dt,f_type,min_distance_Med);
      if( verbose > 1 ) cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance: "<<min_distance_Min<<endl;
      if( verbose > 1 ) cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance: "<<min_distance_Med<<endl;
    }///End Data sample loop
    
    
    t2.Stop();
    delete fData;
    return mdists;
  };
  
  ///Calls analysis through TProcPool
  TProcPool workers(N_Cores);
  auto f_hist = (TH2D*)workers.ProcTree(MCs, workItem);  
  f_hist->GetXaxis()->SetTitle("Data Events");
  for(int mcn=0; mcn<int(MC_Names.size()); mcn++)
    f_hist->GetYaxis()->SetBinLabel(mcn+1,MC_Names[mcn]);
  

  ///_______________________ Compute discriminant from MDMCED _____________________________________________________
  cout<<":: [Distance Computing Time]: "; t.Stop(); t.Print(); t.Continue();
  cout<<"\n::::::::::::::::::::::::::::::::[ Computing discriminant ]::::::::::::::::::::::::::::::::::::"<<endl;
  ///--------------------------------------------------------------------------------------------------------------
  Int_t Event;
  Double_t G_PsbD_MinDist;
  vector<Double_t> PsbD_MinDist;
  TTree *tree = new TTree("FastME","Fast Matrix Element Analysis Results");
  tree->SetDirectory(0);
  tree->Branch("Event",&Event,"Event/I");
  tree->Branch("G_PsbD_MinDist",&G_PsbD_MinDist,"G_PsbD_MinDist/D");
  tree->Branch("PsbD_MinDist","vector<Double_t> PsbD_MinDist",&PsbD_MinDist);
  for(Int_t data=0; data<nData; data++){
    Event = data;
    if( verbose != 0 && data%(nData/10) == 0)
      cout<< Form(":: [Remaining]:   %i Events",nData-data) <<endl;

    Double_t min_dr_sig = f_hist->GetBinContent(data+1,1);
    Double_t min_dr_bkg = 1.E15;
    PsbD_MinDist.clear();
    ///Finds closet MC
    for(Int_t mcs=1; mcs<N_MCT; mcs++){
      PsbD_MinDist.push_back( PsbD(min_dr_sig, f_hist->GetBinContent(data+1,mcs+1)) );
      if( f_hist->GetBinContent(data+1,mcs+1) < min_dr_bkg )
	min_dr_bkg = f_hist->GetBinContent(data+1,mcs+1);
    }

    G_PsbD_MinDist = PsbD(min_dr_sig, min_dr_bkg);
    if( verbose > 1 )
      cout<< Form("GSigMin:   %f\t\tGBkgMin:   %f\t\tGPsbDMinDist:   %f", min_dr_sig, min_dr_bkg, G_PsbD_MinDist) << endl;
    
    
    tree->Fill();
  }
  
  ///________________________________ Stoping timming ________________________________________________________
  cout<<"\n::::::::::::::::::::::::::::::::::::[ Process Finished ]::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: [Analysis Total Time]: "; t.Stop(); t.Print();
  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"<<endl;
  ///---------------------------------------------------------------------------------------------------------

  ///Saving FastME results
  gSystem->Exec("mkdir -p "+out_path);
  TFile *tmp = TFile::Open(out_path+"/"+out_name+".root","recreate");
  f_hist->Write();
  tree->Write();
  tmp->Close();
  
  return 0;
}
///********************************************************************************************************************
