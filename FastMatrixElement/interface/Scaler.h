///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::[ Scaler - Compute variables escale factors ]:::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::[ Code Designer: Miqueias M. de Almeida ]:::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#ifndef Scaler_h
#define Scaler_h


#include "FastMatrixElement/FastMatrixElement/interface/FmeDefinitions.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <cmath>
#include <iomanip>

#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TStopwatch.h>
#include <TCanvas.h>

///Headers to TProcPool
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>



///Compute scale factors to normalize the deltas
void FindScaleFactors(std::vector<std::string> vMCs, TString TTreeName, TString PtBranch, TString EtaBranch,
                      TString PhiBranch, Double_t *ScaledPt, Double_t *ScaledEta, Double_t *ScaledPhi){

  gROOT->SetBatch();
  TCanvas *temp = new TCanvas();
  Double_t pt_sum = 0, eta_sum = 0, phi_sum = 0;
  Int_t total = vMCs.size();
  for(Int_t isample = 0; isample < total; isample++){

    TFile *ftmp = TFile::Open( (TString)vMCs[isample] );
    TTree *ttmp = (TTree*)ftmp->Get( TTreeName );

    if( *ScaledPt < 0){
      TString draw_pt = PtBranch + " >> stackpt";
      ttmp->Draw(draw_pt);
      TH1D *stackpt = (TH1D*)gDirectory->Get("stackpt");
      if( *ScaledPt == -1 ) pt_sum += stackpt->GetMean();
      if( *ScaledPt == -2 ) pt_sum += stackpt->GetBinCenter( stackpt->GetMaximumBin() );
    }
    if( *ScaledEta < 0){
      TString draw_eta = EtaBranch+" >> stacketa";
      ttmp->Draw(draw_eta);
      TH1D *stacketa = (TH1D*)gDirectory->Get("stacketa");
      if( *ScaledEta == -1 ) eta_sum += fabs( stacketa->GetMean() );//Should be used only in case you are in region shifted from 0!
      if( *ScaledEta == -2 ) eta_sum += fabs( stacketa->GetBinCenter(stacketa->GetMinimum()) );
    }
    if( *ScaledPhi < 0){
      TString draw_phi = PhiBranch+" >> stackphi";
      ttmp->Draw(draw_phi);
      TH1D *stackphi = (TH1D*)gDirectory->Get("stackphi");
      if( *ScaledPhi == -1 ) phi_sum += fabs( stackphi->GetMean() );//Should be used only in case you are in region shifted from 0!
      if( *ScaledPhi == -2 ) phi_sum += fabs( stackphi->GetBinCenter(stackphi->GetMinimum()) );
    }

    ftmp->Close();
  }

  if( *ScaledPt  < 0 )  *ScaledPt  = total/pt_sum;
  if( *ScaledEta < 0 )  *ScaledEta = total/eta_sum;
  if( *ScaledPhi < 0 )  *ScaledPhi = total/phi_sum;

  temp->Close();

  return;
}



#endif
