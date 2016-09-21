///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::[ Discriminant Module - Analyze Cartographer's Output ]::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Arbiter_h
#define Arbiter_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
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




///================ Discriminant Based on Distance ===================
Double_t GetPsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}
///===================================================================




///_______________________ Compute discriminant from MDMCED _________________________________________________
void Arbiter(FmeSetup Setup){
  

  Int_t nDtFiles	= Setup.vDatas.size(); //Number of data files
  Int_t nMcFiles	= Setup.vMCs.size();   //Number of MC files
  Int_t verbose		= Setup.Verbose;  

  TStopwatch t2;
  t2.Start();
  std::cout<<ansi_blue<<":::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Computing Discriminant"<<ansi_blue<<" ]::::::::::::::::::::::::::::::::::::"<<ansi_reset<<std::endl;

  ///Set the input tree
  Int_t DtFile;
  std::vector<int> *McFile=0;
  std::vector<double> *Mdist=0;
  TFile *fmeFile = TFile::Open(Setup.FmeFile,"update");
  TTree *mtree = (TTree*)fmeFile->Get("CartographerResults");
  mtree->SetBranchAddress("DataFile",&DtFile);
  mtree->SetBranchAddress("MinDistance",&Mdist);
  mtree->SetBranchAddress("McFile",&McFile);
  Int_t fentries = mtree->GetEntries();

  Int_t fDtFile;
  Double_t Global_PsbDist;
  TTree *ftree = new TTree("Discriminant","Discriminant from FastME");
  ftree->SetDirectory(0);
  ftree->Branch("DataFile",&fDtFile,"DataFile/I");
  ftree->Branch("Global_PsbDist",&Global_PsbDist,"Global_PsbDist/D");

  ///Find the tree sectors (sections of tree for each MC)
  ///The data repeats for every MC type inserted in the analysis
  if(fentries % nMcFiles != 0){
    std::cout<<ansi_red<<"[Error]"<<ansi_reset<<" Something went wrong... non-integer tree sectors!!"<<std::endl;
    std::cout<<"Entries on internal tree: "<<fentries<<std::endl;
    std::cout<<"Number of MCs: "<<nMcFiles<<std::endl;
    std::cout<<"(fentries%nMcFiles): "<< fentries % nMcFiles <<std::endl;
    throw std::exception();
  }
  Int_t TreeSectors[nMcFiles];
  for(Int_t ic=0; ic<(Int_t)nMcFiles; ic++) TreeSectors[ic] = ic*nDtFiles;


  ///To store the local distributions
  TH2D *Local_SigPsbDist = new TH2D("Local_SigPsbDist","Discriminant from Sig to each Bkg",100,-0.05,1.05,nMcFiles,0,nMcFiles);
  Local_SigPsbDist->GetXaxis()->SetTitle("Discriminant");
  Local_SigPsbDist->GetYaxis()->SetTitle("MC file");

  TH2D *Local_BkgPsbDist = new TH2D("Local_BkgPsbDist","Discriminant from Bkg to each Sig",100,-0.05,1.05,nMcFiles,0,nMcFiles);
  Local_BkgPsbDist->GetXaxis()->SetTitle("Discriminant");
  Local_BkgPsbDist->GetYaxis()->SetTitle("MC file");

  TH2D *MinDistances = new TH2D("MinDistances","Minimum Distances",nDtFiles,0,nDtFiles,nMcFiles,0,nMcFiles);
  MinDistances->GetXaxis()->SetTitle("Input file");
  MinDistances->GetYaxis()->SetTitle("MC file");
  
  
  ///Getting results from analysis
  ///Each tree row has all events from a data file
  std::cout<<":: Detected "<<nDtFiles<<" data files..."<<std::endl;
  for(Int_t ifile = 0; ifile < nDtFiles; ifile++){
    ///Goes over the events from current file
    mtree->GetEntry(ifile);    
    const Int_t nevents = Mdist->size();

    ///Detect if signal or background
    bool is_dsig = false;
    bool is_dbkg = false;
    for(Int_t idsig = 0; idsig < (Int_t)Setup.SigData.size(); idsig++){
      if( Setup.SigData.at(idsig) == DtFile ) is_dsig = true;
      //std::cout<<"DtFile/sigData: "<<DtFile<<"/"<<Setup.SigData.at(idsig)<<std::endl;
    }
    for(Int_t idbkg = 0; idbkg < (Int_t)Setup.BkgData.size(); idbkg++){
      if( Setup.BkgData.at(idbkg) == DtFile ) is_dbkg = true;
      //std::cout<<"DtFile/bkgData: "<<DtFile<<"/"<<Setup.BkgData.at(idbkg)<<std::endl;
    }

         if(is_dsig == true && is_dbkg == false) fDtFile = 0;
    else if(is_dsig == false && is_dbkg == true) fDtFile = 1;
    else{
      std::cout<<"[NOTE: Skipping data from file "<<ifile<<"!]"<<std::endl;
      continue;
    }
    


    std::vector<Int_t> sig_mc_file, bkg_mc_file;
    std::vector<Double_t> local_min_dr_sig, local_min_dr_bkg;
    for(Int_t ievent = 0; ievent < nevents; ievent++){
      if(nevents-ievent  == 1){
        std::cout<<":: ["<<ansi_violet<<"Remaning DataFile/Events"<<ansi_reset<<"]:  "<<Form("%i/%i/ %.3fseg",ifile,nevents-ievent,t2.RealTime())<<std::endl;
        t2.Continue();
      }
    

      //Reseting the storage
      sig_mc_file.clear();
      bkg_mc_file.clear();
      local_min_dr_sig.clear();
      local_min_dr_bkg.clear();


      Double_t global_min_dr_sig = 1.e15, global_min_dr_bkg = 1.e15;
      ///Looping over the MC sectors
      for(Int_t ic = 0; ic < (Int_t)nMcFiles; ic++){
        mtree->GetEntry(ifile + TreeSectors[ic]); //TreeSectors aligns the results from different MCs

	if(verbose > 1)
          std::cout<<"Loading entry "<<ifile + TreeSectors[ic]<<"\tDataFile/Event/MCCat "<<ifile<<"/"<<ievent<<"/"<<(*McFile).at(ievent)<<std::endl;
      
      
	///Detect if signal or background
	bool is_sig = false;
	bool is_bkg = false;
        for(Int_t isig = 0; isig < (Int_t)Setup.SigMC.size(); isig++)
	  if( Setup.SigMC.at(isig) == (*McFile).at(ievent) ) is_sig = true;
        for(Int_t ibkg = 0; ibkg < (Int_t)Setup.BkgMC.size(); ibkg++)
          if( Setup.BkgMC.at(ibkg) == (*McFile).at(ievent) ) is_bkg = true;


	///Finds closest MC Signal
        if( is_sig == true ){
          sig_mc_file.push_back( (*McFile).at(ievent) );
          local_min_dr_sig.push_back( (*Mdist).at(ievent) );
   	  if( (*Mdist).at(ievent) < global_min_dr_sig ){
	    global_min_dr_sig = (*Mdist).at(ievent);
	  }
        }

        ///Finds closest MC Background
        if( is_bkg == true ){
          bkg_mc_file.push_back( (*McFile).at(ievent) );
          local_min_dr_bkg.push_back( (*Mdist).at(ievent) );
          if( (*Mdist).at(ievent) < global_min_dr_bkg ){
            global_min_dr_bkg = (*Mdist).at(ievent);
	  }
        }

        
        MinDistances->Fill(DtFile, (*McFile).at(ievent), (*Mdist).at(ievent)/float(nevents));
      }//Ending the full verification for a event

      Global_PsbDist = GetPsbD(global_min_dr_sig, global_min_dr_bkg);
      //if( Global_PsbDist == 0) std::cout<<"@@@@@@@@@ ifile/ievent/min_dr_sig/min_dr_bkg: "<<ifile<<"/"<<ievent<<"/"<<global_min_dr_sig<<"/"<<global_min_dr_bkg<<std::endl;

      if(fDtFile == 0){
	for(Int_t ist = 0; ist < (Int_t)bkg_mc_file.size(); ist++){
	  Local_SigPsbDist->Fill( GetPsbD(global_min_dr_sig, local_min_dr_bkg.at(ist)), bkg_mc_file.at(ist), 1 );
	}
      }
      if(fDtFile == 1){
        for(Int_t ist = 0; ist < (Int_t)sig_mc_file.size(); ist++){
          Local_BkgPsbDist->Fill( GetPsbD(local_min_dr_sig.at(ist), global_min_dr_bkg), sig_mc_file.at(ist), 1 );
        }
      }

      if( verbose > 1 )
        std::cout<< Form("GSigMin:   %f\t\tGBkgMin:   %f\t\tGPsbDMinDist:   %f", global_min_dr_sig, global_min_dr_bkg, Global_PsbDist) << std::endl;

    
      ftree->Fill();
    }//End of loop over events
  }//End of loop over data files
  
  ///________________________________ Stoping timming ________________________________________________________
  std::cout<<ansi_blue<<std::endl;
  std::cout<<"::::::::::::::::::::::::::::::::::::[ "<<ansi_cyan<<"Process Finished"<<ansi_blue<<" ]:::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<":: ["<<ansi_cyan<<"Computing Total Time"<<ansi_blue<<"]: "<<Form("%.3f seg", t2.RealTime())<<std::endl;
  std::cout<<":: ["<<ansi_cyan<<"Producing plot results..."<<ansi_blue<<"]"<<std::endl;
  std::cout<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
  std::cout<<ansi_reset<<std::endl;
  ///--------------------------------------------------------------------------------------------------------



  std::cout<<"Define the path: ";
  TString dir_name;
  std::cin >> dir_name;
  fmeFile->mkdir(dir_name);
  fmeFile->cd(dir_name);

  ///Producing some plots
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

  //----------------------------------------------------------------
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


  ///------------ Graph for significancy of cut ------------------------
  TGraph *sigma = new TGraph();


  ///------------- Plot ROC curve for all MCs -----------------------------  
  Int_t jpoint = 0;
  Double_t cutoff, integral=0, max_sigma = 0, best_cut;
  const int discret = 1000;
  float TPR[discret], FPR[discret], TP, FP, TN, FN;
  std::cout<<"Performing results..."<<std::endl;
  for(int j=0; j<discret; j++){
    cutoff = j/float(discret);
    TP = FP = TN = FN = 0;
    for(Int_t i=0; i<(Int_t)ftree->GetEntries(); i++){
      ftree->GetEntry(i);
	  
      
	if(fDtFile == 0){
	  if( Global_PsbDist > cutoff ) TP++;
	  if( Global_PsbDist < cutoff ) FN++;
	  if( j == 0)
	    hsig->Fill(Global_PsbDist);
	}
      


	if(fDtFile == 1){
	  if( Global_PsbDist > cutoff ) FP++;
	  if( Global_PsbDist < cutoff ) TN++;
	  if( j == 0)
	    hbkg->Fill(Global_PsbDist);
	}
      

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


  //Does not show on fly
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("Discriminant_distributions","",0,0,800,800);
  ///----------- Discriminants plot -------------------
  c1->cd();
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

  c1->cd();
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
  c1->Write();
  c1->Close();



  TCanvas *c2 = new TCanvas("ROC Curve","",0,0,800,800);  
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
  
  //c1->cd(4);
  roc->Draw("AL");
  l3->Draw();
  l50->Draw();
  l80->Draw();
  l90->Draw();
  l95->Draw();
  l100->Draw();
  c2->Write();
  c2->Close();
  ///----------------------------------------------------------------------


  Local_SigPsbDist->Write();
  Local_BkgPsbDist->Write(); 
  MinDistances->Write();
  ftree->Write();
  fmeFile->Close();


  ///Clean memory from used pointers
  delete McFile;
  delete Mdist;

  return;

}


#endif
