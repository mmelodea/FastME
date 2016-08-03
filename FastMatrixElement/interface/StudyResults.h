///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::[ To do Some Studies From the FastME Results]:::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#ifndef StudyResults_h
#define StudyResults_h

#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <cassert>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPaveText.h>
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

  
  ///-------------- Define files to be used -----------------------------------
  int sig_id = 99, bkg_id = 99;
  std::vector<int> sig_files, bkg_files;
  std::cout<<"Which samples from config file should be signal? (insert -1 to finish the list)"<<std::endl;
  while(sig_id != -1){
    std::cin >> sig_id;
    if(sig_id >= 0)
      sig_files.push_back(sig_id);
  }
  std::cout<<"Which samples from config file should be background? (insert -1 to finish the list)"<<std::endl;
  while(bkg_id != -1){
    std::cin >> bkg_id;
    if(bkg_id >= 0)
      bkg_files.push_back(bkg_id);
  }
  
  //list files
  if(sig_files.size() == 0 || bkg_files.size() == 0){
    std::cout<<"Inputs not properly inserted!"<<std::endl;
    throw std::exception();
  }
  std::cout<<"Listing the files to be used: "<<std::endl;
  std::cout<<"For Signal:";
  for(int ifs=0; ifs<(int)sig_files.size(); ifs++)
    std::cout<<"  "<<sig_files.at(ifs);
  std::cout<<"\nFor Background:";
  for(int ifb=0; ifb<(int)bkg_files.size(); ifb++)
    std::cout<<"  "<<bkg_files.at(ifb);
  std::cout<<std::endl;
  
  
  ///-------------- Plot discriminant from singal and background --------------
  TFile *FmeInFile = TFile::Open(UserSetup.FmeFile);
  
  Int_t DtFile;
  std::vector<Int_t> *McIndex=0, *McCat=0, *DtObjFlag=0;
  Double_t Global_PsbDist;
  std::vector<Double_t> *MinDist=0, *Local_PsbDist=0;
  TTree *FmeTree = (TTree*)FmeInFile->Get("FastME");
  FmeTree->SetBranchAddress("DataFile",&DtFile);
  FmeTree->SetBranchAddress("PairedMC",&McIndex);
  FmeTree->SetBranchAddress("PairedMCType",&McCat);
  FmeTree->SetBranchAddress("DataObjFlag",&DtObjFlag);
  FmeTree->SetBranchAddress("MinDistance",&MinDist);
  FmeTree->SetBranchAddress("Global_PsbDist",&Global_PsbDist);
  FmeTree->SetBranchAddress("Local_PsbDist",&Local_PsbDist);
  
  TH1D *hsig = new TH1D("hsig","Discriminant Based on Distance",100,-0.05,1.05);
  hsig->SetLineColor(9);
  hsig->SetFillColor(9);
  hsig->SetFillStyle(3001);
  hsig->GetXaxis()->SetTitle("P_{SB}(Distance)");
  hsig->GetYaxis()->SetTitle("Events/0.01");

  TH1D *hbkg = new TH1D("hbkg","Discriminant Based on Distance",100,-0.05,1.05);
  hbkg->SetLineColor(2);
  hbkg->SetFillColor(2);
  hbkg->SetFillStyle(3001);
  hbkg->GetXaxis()->SetTitle("P_{SB}(Distance)");
  hbkg->GetYaxis()->SetTitle("Events/0.01");
    
  TLegend *leg = new TLegend(0.57,0.77,0.87,0.87);
  leg->AddEntry(hsig,"Signal","f");
  leg->AddEntry(hbkg,"Background","f");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  ///----------------------------------------------------------------------  


  ///-------------- Local Discriminant ------------------------------------
  Int_t Nbins = ((Int_t)sig_files.size() + (Int_t)bkg_files.size());
  Int_t Nfiles = Nbins;
  TH2D *LocalDisc = new TH2D("LocalDisc","Global Discriminant by Data File",100,0,1,Nbins,0,Nfiles);
  LocalDisc->GetXaxis()->SetTitle("Discriminant Value");
  LocalDisc->GetYaxis()->SetTitle("Files");
  
  
  ///------------ Plot DR min correlation ---------------------------------
  TH2D *SigDrMins = new TH2D("SigDrMins","Minimum Distances",100,0,10,100,0,10);
  SigDrMins->GetXaxis()->SetTitle("DR_{min} to Signal");
  SigDrMins->GetYaxis()->SetTitle("DR_{min} Background");
  SigDrMins->SetMarkerColor(9);
  
  TH2D *BkgDrMins = new TH2D("BkgDrMins","Minimum Distances",100,0,10,100,0,10);
  BkgDrMins->GetXaxis()->SetTitle("DR_{min} to Signal");
  BkgDrMins->GetYaxis()->SetTitle("DR_{min} to Background");
  BkgDrMins->SetMarkerColor(2);

  ///------------ Graph for significancy of cut ------------------------
  TGraph *sigma = new TGraph();


  std::vector<Int_t> SSigIndexes, BSigIndexes, SBkgIndexes, BBkgIndexes;
  Int_t max_ssig_index = 0, max_sbkg_index = 0, max_bsig_index = 0, max_bbkg_index = 0;
  ///------------- Plot ROC curve for all MCs -----------------------------  
  Int_t jpoint = 0;
  Double_t cutoff, integral=0, max_sigma = 0, best_cut;
  const int discret = 1000;
  float TPR[discret], FPR[discret], TP, FP, TN, FN;
  std::cout<<"Performing results..."<<std::endl;
  for(int j=0; j<discret; j++){
    cutoff = j/float(discret);
    TP = FP = TN = FN = 0;
    for(Int_t i=0; i<(Int_t)FmeTree->GetEntries(); i++){
      FmeTree->GetEntry(i);
	  
      for(int isg=0; isg<(int)sig_files.size(); isg++){
	if(DtFile == sig_files.at(isg)){
	  if( Global_PsbDist > cutoff ) TP++;
	  if( Global_PsbDist < cutoff ) FN++;
	  if( j == 0){
	    hsig->Fill(Global_PsbDist);
	    for(int ibkg=1; ibkg<(int)McIndex->size(); ibkg++){
	      if((*McIndex).at(0) > max_ssig_index) max_ssig_index = (*McIndex).at(0);
              if((*McIndex).at(ibkg) > max_sbkg_index) max_sbkg_index = (*McIndex).at(ibkg);
	      SSigIndexes.push_back((*McIndex).at(0));
              SBkgIndexes.push_back((*McIndex).at(ibkg));
	      SigDrMins->Fill((*MinDist).at(0),(*MinDist).at(ibkg),1);
	    }
	  }
	}
      }

      for(int ibg=0; ibg<(int)bkg_files.size(); ibg++){
	if(DtFile == bkg_files.at(ibg)){
	  if( Global_PsbDist > cutoff ) FP++;
	  if( Global_PsbDist < cutoff ) TN++;
	  if( j == 0){
	    hbkg->Fill(Global_PsbDist);
	    for(int ibkg=1; ibkg<(int)McIndex->size(); ibkg++){
              if((*McIndex).at(0) > max_bsig_index) max_bsig_index = (*McIndex).at(0);
              if((*McIndex).at(ibkg) > max_bbkg_index) max_bbkg_index = (*McIndex).at(ibkg);
              BSigIndexes.push_back((*McIndex).at(0));
              BBkgIndexes.push_back((*McIndex).at(ibkg));
	      BkgDrMins->Fill((*MinDist).at(0),(*MinDist).at(ibkg),1);
	    }
	  }
	}
      }

      if(j == 0)
	LocalDisc->Fill(Global_PsbDist,DtFile,1);

    }//End loop over events
    TPR[j] = TP/float(TP + FN);
    FPR[j] = FP/float(FP + TN);
    if(j>0)
     integral += fabs(FPR[j-1]-FPR[j])*TPR[j];

    ///Computing significancy of best cut
    if(TP+FP != 0){
      sigma->SetPoint(jpoint, cutoff, sqrt(TPR[j]+FPR[j])-sqrt(FPR[j]));
      jpoint++;
      if( (sqrt(TPR[j]+FPR[j])-sqrt(FPR[j])) > max_sigma ){
        max_sigma = sqrt(TPR[j]+FPR[j])-sqrt(FPR[j]);
	best_cut = cutoff;
	//std::cout<<"MaxS: "<<max_sigma<<"\tCut: "<<cutoff<<std::endl;
      }
    }

  }//End loop over cuts
  

  ///------------ Plot frequency of paired MC events ----------------------
  TH1I *SSigPairedfrequency1D = new TH1I("SSigPairedfrequency1D","",max_ssig_index,0,max_ssig_index+1);
  TH1I *SBkgPairedfrequency1D = new TH1I("SBkgPairedfrequency1D","",max_sbkg_index,0,max_sbkg_index+1);
  //TH2I *SigPairedfrequency = new TH2I("SigPairedfrequency","Signal",10,0,max_ssig_index+1,10,0,max_sbkg_index+1);
  //SigPairedfrequency->GetXaxis()->SetTitle("Paired to Signal");
  //SigPairedfrequency->GetYaxis()->SetTitle("Paired to Background");
  assert(SSigIndexes.size() == SBkgIndexes.size());
  for(Int_t sdx=0; sdx<(Int_t)SSigIndexes.size(); sdx++){
    //SigPairedfrequency->Fill(SSigIndexes.at(sdx),SBkgIndexes.at(sdx),1);
    SSigPairedfrequency1D->Fill(SSigIndexes.at(sdx));
    SBkgPairedfrequency1D->Fill(SBkgIndexes.at(sdx));
  }

  TH1I *BSigPairedfrequency1D = new TH1I("BSigPairedfrequency1D","",max_bsig_index,0,max_bsig_index+1);
  TH1I *BBkgPairedfrequency1D = new TH1I("BBkgPairedfrequency1D","",max_bbkg_index,0,max_bbkg_index+1);
  //TH2I *BkgPairedfrequency = new TH2I("BkgPairedfrequency","Background",10,0,max_bsig_index+1,10,0,max_bbkg_index+1);
  //BkgPairedfrequency->GetXaxis()->SetTitle("Paired to Signal");
  //BkgPairedfrequency->GetYaxis()->SetTitle("Paired to Background");
  assert(BSigIndexes.size() == BBkgIndexes.size());
  for(Int_t bdx=0; bdx<(Int_t)BSigIndexes.size(); bdx++){
    //BkgPairedfrequency->Fill(BSigIndexes.at(bdx),BBkgIndexes.at(bdx),1);
    BSigPairedfrequency1D->Fill(BSigIndexes.at(bdx));
    BBkgPairedfrequency1D->Fill(BBkgIndexes.at(bdx));
  }



  //Does not show on fly
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TString plot_name;
  TCanvas *c1 = new TCanvas("c1","",0,0,1280,1024);
  c1->Divide(3,2);  

  
  ///----------- Discriminants plot -------------------
  c1->cd(1);
  TPad* pad1 = new TPad("pad1","1",0,0.4,1,0.95); pad1->Draw(); pad1->cd();
  pad1->SetBottomMargin(0.);
  if(hsig->GetMaximum() > hbkg->GetMaximum()){
    hsig->Draw();
    hbkg->Draw("same");
  }
  else{
    hbkg->Draw();
    hsig->Draw("same");
  }
  leg->Draw();

  c1->cd(1);
  TPad* pad2 = new TPad("pad2","2",0,0.05,1,0.396); pad2->Draw(); pad2->cd();
  pad2->SetTopMargin(0.);
  pad2 -> SetGridx(true);
  pad2 -> SetGridy(true);
  sigma->Draw("AP");
  sigma->GetXaxis()->SetLimits(-0.05,1.05);
  sigma->GetXaxis()->SetTitle("Cut Value");
  sigma->GetXaxis()->SetTitleSize(0.1);
  sigma->GetXaxis()->SetLabelSize(0.09);
  sigma->GetYaxis()->SetTitle("Significance");
  sigma->GetYaxis()->SetTitleSize(0.1);
  sigma->GetYaxis()->SetTitleOffset(0.2);
  sigma->GetYaxis()->SetLabelSize(0.);
  sigma->SetLineColor(kViolet);
  sigma->SetMarkerColor(kViolet);

  TPaveText sigInfo(0.7,0.8,0.9,0.9,"NDC");
  sigInfo.SetFillStyle(0);
  sigInfo.SetBorderSize(0);
  sigInfo.AddText(Form("Best: %.3f",best_cut));
  sigInfo.Draw();
  gPad->Update();  
  
  ///--------------- ROC plot ---------------------------------------
  std::cout<<Form("Area under ROC curve = %.3f",integral)<<std::endl;
  std::cout<<Form("Best significance at cut = %.3f",best_cut)<<std::endl;
  
  TGraph *roc = new TGraph(discret,FPR,TPR);
  roc->SetTitle(Form("ROC Plot - Area under curve = %.3f",integral));
  roc->SetLineColor(kOrange);
  roc->SetLineWidth(2);
  roc->GetXaxis()->SetTitle("False Positive Rate");
  roc->GetXaxis()->SetRangeUser(0,1.);
  roc->GetYaxis()->SetTitle("True Positive Rate");
  roc->GetYaxis()->SetRangeUser(0,1.1);
  
  TLine *l3 = new TLine(0,0,1,1);
  l3->SetLineColor(kGreen-6);
  l3->SetLineStyle(2);

  TLine *l50  = new TLine(0,0.5,1,0.5);		l50->SetLineStyle(2);
  TLine *l80  = new TLine(0,0.8,1,0.8);         l80->SetLineStyle(2);
  TLine *l90  = new TLine(0,0.9,1,0.9);         l90->SetLineStyle(2);
  TLine *l95  = new TLine(0,0.95,1,0.95);       l95->SetLineStyle(2);
  TLine *l100 = new TLine(0,1.,1,1.);           l100->SetLineStyle(2);
  
  c1->cd(4);
  roc->Draw("AL");
  l3->Draw();
  l50->Draw();
  l80->Draw();
  l90->Draw();
  l95->Draw();
  l100->Draw();
  ///----------------------------------------------------------------------
  
  
  c1->cd(2);
  TPad* pad7 = new TPad("pad7","7",0,0.55,1,0.95); pad7->Draw(); pad7->cd();
  pad7->SetBottomMargin(0.);
  pad7->SetGridx(true);
  pad7->SetGridy(true);
  SigDrMins->Draw("Colz");
  c1->cd(2);
  TPad* pad8 = new TPad("pad8","8",0,0.05,1,0.5); pad8->Draw(); pad8->cd();
  pad8->SetTopMargin(0.);
  pad8->SetGridx(true);
  pad8->SetGridy(true);
  BkgDrMins->Draw("Colz");


  c1->cd(5);
  LocalDisc->Draw("Colz");


  c1->cd(3);
  TPad* pad3 = new TPad("pad3","3",0,0.55,1,0.95); pad3->Draw(); pad3->cd();
  pad3->SetBottomMargin(0.);
  SSigPairedfrequency1D->SetFillColor(kBlue);
  SSigPairedfrequency1D->SetLineColor(kBlue);
  SSigPairedfrequency1D->SetFillStyle(3002);
  SSigPairedfrequency1D->GetYaxis()->SetTitle("Signal Frequency");
  SBkgPairedfrequency1D->SetFillColor(kRed);
  SBkgPairedfrequency1D->SetLineColor(kRed);
  SBkgPairedfrequency1D->SetFillStyle(3002);
  SBkgPairedfrequency1D->GetYaxis()->SetTitle("Background Frequency");
  SBkgPairedfrequency1D->GetXaxis()->SetTitle("Event Index");
  SSigPairedfrequency1D->Draw();
  c1->cd(3);
  TPad* pad4 = new TPad("pad4","4",0,0.05,1,0.5); pad4->Draw(); pad4->cd();
  pad4->SetTopMargin(0.);
  SBkgPairedfrequency1D->Draw();


  
  c1->cd(6);
  TPad* pad5 = new TPad("pad5","5",0,0.55,1,0.95); pad5->Draw(); pad5->cd();
  pad5->SetBottomMargin(0.);
  BSigPairedfrequency1D->SetFillColor(kBlue);
  BSigPairedfrequency1D->SetLineColor(kBlue);
  BSigPairedfrequency1D->SetFillStyle(3002);
  BSigPairedfrequency1D->GetYaxis()->SetTitle("Signal Frequency");
  BBkgPairedfrequency1D->SetFillColor(kRed);
  BBkgPairedfrequency1D->SetLineColor(kRed);
  BBkgPairedfrequency1D->SetFillStyle(3002);
  BBkgPairedfrequency1D->GetYaxis()->SetTitle("Background Frequency");
  BBkgPairedfrequency1D->GetXaxis()->SetTitle("Event Index");
  BSigPairedfrequency1D->Draw();
  c1->cd(6);
  TPad* pad6 = new TPad("pad6","6",0,0.05,1,0.5); pad6->Draw(); pad6->cd();
  pad6->SetTopMargin(0.);
  BBkgPairedfrequency1D->Draw();


  c1->Update();
  std::cout<<"Enter name to save plot: "; 
  std::cin >> plot_name;
  c1->Print(UserSetup.OutPath+"/"+plot_name+".png");



  delete McIndex;
  delete McCat;
  delete DtObjFlag;
  delete MinDist;
  delete Local_PsbDist;
  
  
  return;
}


#endif
