///To see in how are separated the particles in the (pT,eta,phi) phase space
void display_particles(){
  gROOT->SetBatch(kTRUE);
  
  TFile *fmeResults = TFile::Open("FME_results/fme_results_data.root");
  Double_t RParticlePt[4], RParticleEta[4];  
  TH2D *mdists = (TH2D*)fmeResults->Get("mdists");
  TH2D *indices = (TH2D*)fmeResults->Get("indices");
  Int_t nentries = indices->GetNbinsX();
  
  TFile *inputs[4];
  inputs[0] = TFile::Open("CMS_Ntuples/FastME_format/4l_data.root");
  inputs[1] = TFile::Open("CMS_Ntuples/FastME_format/4l_SM_Higgs126_2.root");
  inputs[2] = TFile::Open("CMS_Ntuples/FastME_format/4l_ggZZ_2.root");
  inputs[3] = TFile::Open("CMS_Ntuples/FastME_format/4l_qqZZ_2.root");
  
  Double_t OParticlePt[4], OParticleEta[4], OParticlePhi[4], mindr;
  Int_t evIds[7], color[4] = {1,4,2,6}, marker[4] = {28,24,3,26};
  TFile *distances = new TFile("particles_distances_display.root","recreate");
  //nentries = 1;
  TMultiGraph *full[3];
  TGraph *phs[4];
  TString view[3] = {"eta_pt_view","phi_pt_view","eta_phi_view"};
  cout<<"Getting event ";
  TCanvas *cv = new TCanvas("cv","cv",0,0,500,500);
  TCanvas *cv2 = new TCanvas("joint_view","Joint View",0,400,1000,400);
  cv2->Divide(3,1);
  for(Int_t i=0; i<nentries; i++){
    cout<<i<<", ";
    distances->mkdir(Form("Event %i",i));
    distances->cd(Form("Event %i",i));

    evIds[0] = i;///For data
    
    evIds[1] = indices->GetBinContent(i+1,1);
    evIds[2] = indices->GetBinContent(i+1,2);
    evIds[3] = indices->GetBinContent(i+1,3);
    
    ///look to minimum dr
    if(i==0){
      mindr = (mdists->GetBinContent(i+1,1) < mdists->GetBinContent(i+1,2))? mdists->GetBinContent(i+1,1) : mdists->GetBinContent(i+1,2);
      mindr = (mindr < mdists->GetBinContent(i+1,3))? mindr : mdists->GetBinContent(i+1,3);
    }
    if(i>0){
      mindr = (mindr < mdists->GetBinContent(i+1,1))? mindr : mdists->GetBinContent(i+1,1);
      mindr = (mindr < mdists->GetBinContent(i+1,2))? mindr : mdists->GetBinContent(i+1,2);
      mindr = (mindr < mdists->GetBinContent(i+1,3))? mindr : mdists->GetBinContent(i+1,3);
    }
    
    ///3 views
    full[0] = new TMultiGraph();
    full[1] = new TMultiGraph();
    full[2] = new TMultiGraph();
    for(Int_t s=0; s<4; s++){
      TTree *tree = (TTree*)inputs[s]->Get("Higgs");
      tree->SetBranchAddress("ParticlePt",&OParticlePt);
      tree->SetBranchAddress("ParticleEta",&OParticleEta);
      tree->SetBranchAddress("ParticlePhi",&OParticlePhi);
      tree->GetEntry(evIds[s]);

      for(Int_t g=0; g<4; g++){
	phs[g] = new TGraph();
	phs[g]->SetMarkerStyle(marker[s]);
	phs[g]->SetMarkerSize(1.4);
	phs[g]->SetMarkerColor(color[s]);
      }
      Int_t ipoint = 0;
      for(Int_t p=0; p<4; p++){
	  phs[0]->SetPoint(ipoint,OParticlePt[p],OParticleEta[p]);
	  phs[1]->SetPoint(ipoint,OParticlePt[p],OParticlePhi[p]);
	  phs[2]->SetPoint(ipoint,OParticleEta[p],OParticlePhi[p]);
	  ipoint++;
	}
	full[0]->Add( phs[0], "p" );
	full[1]->Add( phs[1], "p" );
	full[2]->Add( phs[2], "p" );
    }

    cv->cd();
    full[0]->Draw("ap");
    full[0]->SetName(view[0]);
    full[0]->SetTitle("#eta vs. p_{T}[GeV]");
    full[0]->GetXaxis()->SetTitle("p_{T}[GeV]");
    full[0]->GetYaxis()->SetTitle("#eta");
    cv->Update();
    full[0]->Write();
    
    full[1]->Draw("ap");
    full[1]->SetName(view[1]);
    full[1]->SetTitle("#phi vs. p_{T}[GeV]");
    full[1]->GetXaxis()->SetTitle("p_{T}[GeV]");
    full[1]->GetYaxis()->SetTitle("#phi");
    cv->Update();
    full[1]->Write();

    full[2]->Draw("ap");
    full[2]->SetName(view[2]);
    full[2]->SetTitle("#eta vs. #phi");
    full[2]->GetXaxis()->SetTitle("#phi");
    full[2]->GetYaxis()->SetTitle("#eta");
    cv->Update();
    full[2]->Write();
    
    cv2->cd(1);
    full[0]->Draw("ap");
    full[0]->SetTitle("#eta vs. p_{T}[GeV]");
    full[0]->GetXaxis()->SetTitle("p_{T}[GeV]");
    full[0]->GetYaxis()->SetTitle("#eta");
    cv2->Update();
    cv2->cd(2);
    full[1]->Draw("ap");
    full[1]->SetTitle("#phi vs. p_{T}[GeV]");
    full[1]->GetXaxis()->SetTitle("p_{T}[GeV]");
    full[1]->GetYaxis()->SetTitle("#phi");
    cv2->Update();
    cv2->cd(3);
    full[2]->Draw("ap");
    full[2]->SetTitle("#eta vs. #phi");
    full[2]->GetXaxis()->SetTitle("#phi");
    full[2]->GetYaxis()->SetTitle("#eta");
    gPad->Modified();
    cv2->Update();
    cv2->Write();
  }
  distances->Close();
  cout<<endl;
  cout<<"Minimum DR = "<<mindr<<endl;
}//End script
