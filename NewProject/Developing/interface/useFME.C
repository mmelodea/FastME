///::::::::::::::::::::::::: To Use the FME Results :::::::::::::::::::::::
///:::::::::::::::::::::: Author: Miquéias M. Almeida :::::::::::::::::::::

int useFME(){

  TString path[3][2];
//  path[0][0] = "/home/sabayon/GitHub/Ntuples/higgs_data_formated.root";
//  path[0][1] = "Data_FME_Results.root";
  
  path[1][0] = "/home/sabayon/GitHub/Ntuples/VBF_QED4_QCD0_SIG1.root";
  path[1][1] = "Sig_FME_Results.root";
  
  path[2][0] = "/home/sabayon/GitHub/Ntuples/VBF_QED4_QCD0_BKG1.root";
  path[2][1] = "Bkg_FME_Results.root";
  
  TH2D *mc_sig = new TH2D("mc_sig","",70,100,800,50,0,1);
  //TH2D *mc_sig = new TH2D("mc_sig","",54,12,120,50,0,1);
  mc_sig->SetMarkerColor(kBlue);
  mc_sig->SetFillColor(kBlue);
  TH2D *mc_bkg = new TH2D("mc_bkg","",70,100,800,50,0,1);
  //TH2D *mc_bkg = new TH2D("mc_bkg","",54,12,120,50,0,1);
  mc_bkg->SetMarkerColor(kRed);
  mc_bkg->SetFillColor(kRed);
  TH2D *data   = new TH2D("data","",70,100,800,50,0,1);
  //TH2D *data   = new TH2D("data","",54,12,120,50,0,1);
  data->SetMarkerStyle(26);
  
  for(int s=1; s<3; s++){
    Float_t f_mass;
    TFile *file = TFile::Open(path[s][0]);
    TTree *file_Tree = (TTree*)file->Get("Higgs");
    file_Tree->SetDirectory(0);
    //file_Tree->SetBranchAddress("Zon_mass",&f_mass);
    //file_Tree->SetBranchAddress("Zoff_mass",&f_mass);
    file_Tree->SetBranchAddress("H_reco_mass",&f_mass);
  
    double prob_sig_bkg, ws, wb, w;
    TFile *fFME = TFile::Open(path[s][1]);
    TTree *FME_Tree = (TTree*)fFME->Get("FastME_Results");
    FME_Tree->SetDirectory(0);
    int nFME = FME_Tree->GetEntries();
    FME_Tree->SetBranchAddress("P_SB",&prob_sig_bkg);
    FME_Tree->SetBranchAddress("WSig_ToEvent",&ws);
    FME_Tree->SetBranchAddress("WBkg_ToEvent",&wb);
    FME_Tree->SetBranchAddress("Event_Weight",&w);
  
    for(int i=0; i<nFME; i++){
      file_Tree->GetEntry(i);
      FME_Tree->GetEntry(i);
      if(s==0) data->Fill(f_mass,prob_sig_bkg);
      if(s==1) mc_sig->Fill(f_mass,prob_sig_bkg);
      if(s==2) mc_bkg->Fill(f_mass,prob_sig_bkg);
    }
  }
  
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv","",10,10,1500,1300);
  
  mc_bkg->Draw("Colz");
  mc_sig->Draw("Colz,same");
  //data->Draw("P,same");
  
  TLine *l = new TLine(100,0.5,800,0.5);
  l->SetLineColor(kGreen);
  l->SetLineWidth(2);
  l->SetLineStyle(2);
  l->Draw();

  TLegend *leg = new TLegend(0.6,0.4,0.9,0.2);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  //leg->AddEntry(data,"Data","p");
  leg->AddEntry(mc_sig,"pp > H > ZZjj","f");
  leg->AddEntry(mc_bkg,"pp > ZZjj /H","f");  
  leg->AddEntry(l,"Method Cut","l");
  leg->Draw();
  
  return 0; //if well finished
}
