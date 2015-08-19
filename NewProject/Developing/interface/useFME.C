///::::::::::::::::::::::::: To Use the FME Results :::::::::::::::::::::::
///:::::::::::::::::::::: Author: Miquéias M. Almeida :::::::::::::::::::::

int useFME(){

  Float_t f_mass4l;
  TFile *fData = TFile::Open("/home/sabayon/temp/FastME/Descriminador/Higgs/prepare/higgs_data_formated.root");
  TTree *Data_Tree = (TTree*)fData->Get("VBF");
  Data_Tree->SetBranchAddress("reco",&f_mass4l);

  double prob_sig_bkg, ws, wb, w;
  TFile *fFME = TFile::Open("FastME_Results.root");
  TTree *FME_Tree = (TTree*)fFME->Get("FastME_Results");
  int nFME = FME_Tree->GetEntries();
  FME_Tree->SetBranchAddress("P_SB",&prob_sig_bkg);
  FME_Tree->SetBranchAddress("WSig_ToEvent",&ws);
  FME_Tree->SetBranchAddress("WBkg_ToEvent",&wb);
  FME_Tree->SetBranchAddress("Event_Weight",&w);
  
  TH2D *sp4l = new TH2D("sp4l","",70,100,800,50,0,1);
  for(int i=0; i<nFME; i++){
    Data_Tree->GetEntry(i);
    FME_Tree->GetEntry(i);
    sp4l->Fill(f_mass4l,prob_sig_bkg);
  }
  
  sp4l->Draw();
  
  return 0; //if well finished
}