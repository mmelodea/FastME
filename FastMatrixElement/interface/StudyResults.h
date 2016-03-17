///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::[ To do Some Studies From the FastME Results]:::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef StudyResults_h
#define StudyResults_h

#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>


void StudyResults(FmeSetup UserSetup){
  //To get better colors
  Int_t palette[5];
  palette[0] = 1;
  palette[1] = 2;
  palette[2] = 3;
  palette[3] = 4;
  palette[4] = 5;
  palette[5] = 6;
  palette[6] = 7;
  palette[7] = 8;
  palette[8] = 9;
  palette[9] = 10;
  gStyle->SetPalette(10,palette);

  const Int_t Number = 3;
  Double_t Red[Number]    = { 1.00, 0.00, 0.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  
  //Does not show on fly
  gROOT->SetBatch(kTRUE);

  
  ///Plot discriminant from singal and background 
  //TFile *fdata = TFile::Open(UserSetup.OutPath+"/"+UserSetup.OutName+".root");
  TFile *fsig  = TFile::Open(UserSetup.OutPath+"/"+UserSetup.OutName+".root");
  TFile *fbkg  = TFile::Open(UserSetup.OutPath+"/"+UserSetup.OutName+".root");
  
  //TTree *tdata = (TTree)fdata->Get("FastME");
  TTree *tsig  = (TTree*)fsig->Get("FastME");
  TTree *tbkg  = (TTree*)fbkg->Get("FastME");
  
  TH1D *hsig = new TH1D("hsig","hsig",100000,0,1);
  hsig->SetLineColor(9);
  hsig->GetXaxis()->SetTitle("P_{SB}(Distance)");
  hsig->GetYaxis()->SetTitle("Events/0.02 (Normalized)");

  TH1D *hbkg = new TH1D("hbkg","hbkg",100000,0,1);	hbkg->SetLineColor(2);
    
  TCanvas *c1 = new TCanvas("c1","",0,0,600,600);
  tsig->Draw("Global_PsbDist >> hsig");
  tbkg->Draw("Global_PsbDist >> hbkg");
  
  hsig->Draw("9");
  hbkg->Draw("9,same");
  c1->Update();
  c1->Print("Discriminant_Signal_vs_Background.png");
  
/*
  Double_t cutoff, integral=0;
  const int discret = 50;
  float TPR[discret], FPR[discret], TP, FP, TN, FN;
  for(int j=0; j<discret; j++){
    cutoff = j/float(discret);
    TP = FP = TN = FN = 0;
    for(int i=0; i<nevents; i++){
      tsig->GetEntry(i);
      if( GSigPsbDistance > cutoff ) TP++;
      if( GSigPsbDistance < cutoff ) FN++;
      tbkg->GetEntry(i);
      if( GBkgPsbDistance > cutoff ) FP++;
      if( GBkgPsbDistance < cutoff ) TN++;
    }
    TPR[j] = TP/float(TP + FN);
    FPR[j] = FP/float(FP + TN);
    if(j>0)
     integral += fabs(FPR[j-1]-FPR[j])*TPR[j];
    //if(TPR[j]>=0.8) cout<<"TPR: "<<TPR[j]<<"\tFPR: "<<FPR[j]<<"\tCutOff: "<<cutoff<<endl;
  }
  cout<<"Area in the ROC plot: "<<integral<<endl;

  TGraph *roc = new TGraph(discret,FPR,TPR);
  roc->SetTitle("ROC");
  roc->SetMarkerStyle(4);
  roc->SetMarkerSize(0.9);
  roc->SetMarkerColor(kOrange);
  roc->GetXaxis()->SetTitle("False Positive Rate");
  roc->GetXaxis()->SetRangeUser(0,1.);
  roc->GetYaxis()->SetTitle("True Positive Rate");
  roc->GetYaxis()->SetRangeUser(0,1.1);
  
  TLine *l3 = new TLine(0,0,1,1);
  l3->SetLineColor(kBlack);
  l3->SetLineWidth(2);
  l3->SetLineStyle(2);

  TLine *l50  = new TLine(0,0.5,1,0.5);		l50->SetLineStyle(2);
  TLine *l80  = new TLine(0,0.8,1,0.8);         l80->SetLineStyle(2);
  TLine *l90  = new TLine(0,0.9,1,0.9);         l90->SetLineStyle(2);
  TLine *l95  = new TLine(0,0.95,1,0.95);         l95->SetLineStyle(2);
  TLine *l100 = new TLine(0,1.,1,1.);           l100->SetLineStyle(2);
  
  TLegend *leg = new TLegend(0.2,0.75,0.5,0.85);
  leg->AddEntry(sig_psbD,"SM Higgs 126GeV","f");
  leg->AddEntry(bkg_psbD,"ggZZ + qqZZ","f");
  leg->AddEntry(bkg1_psbD,"ggZZ","l");
  leg->AddEntry(bkg2_psbD,"qqZZ","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);  
  
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,1);
  cv->cd(1);
  sig_psbD->DrawNormalized();
  bkg_psbD->DrawNormalized("same");
  bkg1_psbD->DrawNormalized("same");
  bkg2_psbD->DrawNormalized("same");
  leg->Draw();
  pv->Draw();
  pv2->Draw();
  cv->cd(2);
  roc->Draw("AP");
  l3->Draw();
  l50->Draw();
  l80->Draw();
  l90->Draw();
  l95->Draw();
  l100->Draw();
*/
  
}


#endif