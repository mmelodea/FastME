//By Miquéias M. de Almeida 

#include "FastME.h"
#include <time.h>
#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>

#define cut 0.5					//Threshold to separate events (ideal cut 0.5 - MC #Sig and #Bkg equals)

using namespace std;

///****************************************************************************************************///
/// Calculate the distance between events using matrix with columms containing pT, eta, phi and lines  ///
/// containing (l1-l4,j1,j2) the leptons indice. The leptons Data and MC are organized by closure      ///
/// constraining that 2 leptons came from the same vertex (constraining Ressonance)      	       ///
///****************************************************************************************************///
double FME::compute_DR(string type, TTree *Data_tree, TTree *MC_tree){
  switch(type){
    case "4l":
      return FS4l_DR(Float_t *Data,Float_t *MC);
    
    case "4l2j":
      return FS4l2j_DR(Float_t *Data,Float_t *MC);
      
    case "lv2j":
      return FSlv2j_DR(Float_t *Data,Float_t *MC);
    
    default:
      return FS4l_DR(Float_t *Data,Float_t *MC);
  }
}

///Final State 4l only
double FME::FS4l_DR(Float_t *Data, Float_t *MC){
    
    double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
    double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
    int start = 0;
  
    do{
	sum_dr1 = 0;
	sum_dr2 = 0;
      
	fdpt2    = pow(Data[start][0]-MC[start][0],2);//*Data[start][0][1];
	fdeta2   = pow(Data[start][1]-MC[start][1],2);//*Data[start][1][1];
	fdphi2   = pow(Data[start][2]-MC[start][2],2);//*Data[start][2][1];
	sum_dr1 = sqrt(dpt2 + deta2 + dphi2);
      
	sdpt2    = pow(Data[start][0]-MC[start+1][0],2);//*Data[start][0][1];
	sdeta2   = pow(Data[start][1]-MC[start+1][1],2);//*Data[start][1][1];
	sdphi2   = pow(Data[start][2]-MC[start+1][2],2);//*Data[start][2][1];
	sum_dr2 = sqrt(dpt2 + deta2 + dphi2);
      
	if(sum_dr1 < sum_dr2){
	  sum_dpt2  += fdpt2;
	  sum_deta2 += fdeta2;
	  sum_dphi2 += fdphi2;
	
	  sum_dpt2  += pow(Data[start+1][0]-MC[start+1][0],2);
	  sum_deta2 += pow(Data[start+1][1]-MC[start+1][1],2);
	  sum_dphi2 += pow(Data[start+1][2]-MC[start+1][2],2);
	}
      
	if(sum_dr2 < sum_dr1){
	  sum_dpt2  += sdpt2;
	  sum_deta2 += sdeta2;
	  sum_dphi2 += sdphi2;
	
	  sum_dpt2  += pow(Data[start+1][0]-MC[start][0],2);
	  sum_deta2 += pow(Data[start+1][1]-MC[start][1],2);
	  sum_dphi2 += pow(Data[start+1][2]-MC[start][2],2);
	}

	start += 2; 
    }while(start < 4);
  
    event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
    return event_distance;
}

///semi-leptonic case
double FME::FS_lv2j_DR(Float_t *Data, Float_t *MC){
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start=0;

    dpt  = (Data[start][0][0] - MC[start][0][0]);//*pt_res(start,MC[start][0][1],MC[start][1][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//*(1./MC[start][1][1]);
    dphi = (Data[start][2][0] - MC[start][2][0]);//*(1./MC[start][2][1]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    //MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0]);//*(1./MC[start][0][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//*(1./MC[start][1][1]);
    dphi = (Data[start][2][0] - MC[start][2][0]);//*(1./MC[start][2][1]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    //Jet 1
    start = 2;
    dpt  = (Data[start][0][0] - MC[start][0][0]);//*pt_res(start,MC[start][0][1],MC[start][1][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//MC[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0]);//MC[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    //Jet 2
    start = 3;
    dpt  = (Data[start][0][0] - MC[start][0][0]);//*pt_res(start,MC[start][0][1],MC[start][1][1]);
    deta = (Data[start][1][0] - MC[start][1][0]);//MC[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0]);//MC[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///====================================================================================================