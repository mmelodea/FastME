//By Miquéias M. de Almeida

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>

#define l 6					//Number of particles
#define c 4					//Number of variables
#define sigma 1					//Measurement Resolution
#define phspace_radius 100.			//Radius around the event


/*****************************************************************************************************
 * Calculate the distance between events using matrix with columms containing pT, eta, phi and lines *
 * containing l1-l4,j1,j2 the leptons comming from ntuples are organized by closure between them     *
 * (the particle to compare are ordered by the same distance critery)                                *
 *****************************************************************************************************/
Float_t dr(Float_t Data[l][3], Float_t MC[l][c]){
  
  Float_t dpt2, deta2, dphi2, sum_dr1, sum_dr2;
  Float_t sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  //Choose the l/j1-Data most close to l/j1-MC
  do{
      sum_dr1 = 0;
      sum_dr2 = 0;
      
      dpt2    = pow(Data[start][0]-MC[start][0],2);
      deta2   = pow(Data[start][1]-MC[start][1],2);
      dphi2   = pow(Data[start][2]-MC[start][2],2);
      sum_dr1 = sqrt(dpt2 + deta2 + dphi2);
      
      dpt2    = pow(Data[start][0]-MC[start+1][0],2);
      deta2   = pow(Data[start][1]-MC[start+1][1],2);
      dphi2   = pow(Data[start][2]-MC[start+1][2],2);
      sum_dr2 = sqrt(dpt2 + deta2 + dphi2);
      
      if(sum_dr1 < sum_dr2){
	sum_dpt2  += pow(Data[start][0]-MC[start][0],2);
	sum_deta2 += pow(Data[start][1]-MC[start][1],2);
	sum_dphi2 += pow(Data[start][2]-MC[start][2],2);
	
	sum_dpt2  += pow(Data[start+1][0]-MC[start+1][0],2);
	sum_deta2 += pow(Data[start+1][1]-MC[start+1][1],2);
	sum_dphi2 += pow(Data[start+1][2]-MC[start+1][2],2);
      }
      
      if(sum_dr2 < sum_dr1){
	sum_dpt2  += pow(Data[start][0]-MC[start+1][0],2);
	sum_deta2 += pow(Data[start][1]-MC[start+1][1],2);
	sum_dphi2 += pow(Data[start][2]-MC[start+1][2],2);
	
	sum_dpt2  += pow(Data[start+1][0]-MC[start][0],2);
	sum_deta2 += pow(Data[start+1][1]-MC[start][1],2);
	sum_dphi2 += pow(Data[start+1][2]-MC[start][2],2);
      }

    start += 2;
  }while(start < l);
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
//====================================================================================================




void FastME(){
  
  int mc_neighbors, nevents_weighted = 0, chosen_id, neighbors = 0;
  Float_t data_stored[l][3], mc_stored[l][c], mc_chosen[l][3];
  Float_t mc_test, mc_min_distance, min_distance = 1000, mc_max_distance, max_distance = 0, data_weight, dr_sum;
    
  TFile *data_file = TFile::Open("Ntuples/Normal_Data_DR.root");
  TTree *data = (TTree*)data_file->Get("VBF");
  data->SetBranchAddress("RECO_PARTICLE",&data_stored);
  Int_t ndata = (Int_t)data->GetEntries();
  
  TFile *mc_file = TFile::Open("Ntuples/Sig.root");
  TTree *mc = (TTree*)mc_file->Get("VBF");
  mc->SetBranchAddress("RECO_PARTICLE",&mc_stored);
  Int_t nmc = (Int_t)mc->GetEntries();
  
  TTree *results = new TTree("results","Results from FastME");
  results->SetDirectory(0);
  results->Branch("data_weight",&data_weight,"data_weight/F");
  //results->Branch("mc_chosen",&mc_chosen,"mc_chosen[6][3]/F");
    
  //******* Export analysis to txt file *************************************************
  FILE *AnalysisInfo;
  AnalysisInfo = fopen("AnalysisInfo.txt","w");
  fprintf(AnalysisInfo,"%s\n","_________________________________________________________");
  fprintf(AnalysisInfo,"%s\n","              Fast Matrix Element Results                "); 
  fprintf(AnalysisInfo,"%s\n","_________________________________________________________");
  fprintf(AnalysisInfo,"%s\n","---------------------------------------------------------");
  fprintf(AnalysisInfo,"%s\t %i\n","Data Events:",ndata);
  fprintf(AnalysisInfo,"%s\t %i\n","MC   Events:",nmc);
  fprintf(AnalysisInfo,"%s\t %f\n","PHS  RADIUS:",phspace_radius);
  fprintf(AnalysisInfo,"%s\n","_________________________________________________________");
  fprintf(AnalysisInfo,"%s\n","---------------------------------------------------------");
  //*************************************************************************************
  
  printf("\n%s\t %i\n","Data Events: ",ndata);
  printf("%s\t %i\n","MC   Events: ",nmc);
  printf("%s\t %f\n\n","PHS  Radius: ",phspace_radius);
  
  //ndata = 1.E3;
  //nmc = 1.E3;
  for(Int_t i=0; i<ndata; i++){
    data->GetEntry(i);
    if(i % 1000 == 0) 
      cout<<"Remaining Data: "<<ndata-i<<endl;

    //Reset variables
    mc_neighbors    = 0;
    data_weight     = 0;
    mc_min_distance = 1.E9;
    mc_max_distance = 0;
    
    //Choose the most close MC event to current Data event
    dr_sum = 1.E-2;
    for(Int_t j=0; j<nmc; j++){
      mc->GetEntry(j);
      mc_test = dr(data_stored,mc_stored);      
      if(mc_test < phspace_radius){
	mc_neighbors += 1;
	dr_sum += mc_test;
      }
      //if(mc_neighbors > neighbors){ chosen_id = i; neighbors = mc_neighbors; }
      if(mc_test < mc_min_distance) mc_min_distance = mc_test;
      if(mc_test > mc_max_distance) mc_max_distance = mc_test;
    }    
      //Getting the data event weight
      data_weight = mc_neighbors/dr_sum;
      if(mc_neighbors == 0 && dr_sum == 0) data_weight = 0;

      results->Fill();
     
    if(data_weight != 0){
      nevents_weighted += 1;
      fprintf(AnalysisInfo,"%s\t\t %i\n","Data Event:",i);
      fprintf(AnalysisInfo,"%s\t\t %i\n","MC Neighbors:",mc_neighbors);
      fprintf(AnalysisInfo,"%s\t\t\t %f\n","Weight:",data_weight);
      fprintf(AnalysisInfo,"%s\t %f\n\n","Minimum Distance:",mc_min_distance);
    }
      
    //Checks the minimum distance between Data and MC evend found in all process
    if(mc_min_distance < min_distance) min_distance = mc_min_distance;
    if(mc_max_distance > max_distance) max_distance = mc_max_distance;
  }
  
  float fev_weighted = 100*nevents_weighted/float(ndata);
  fprintf(AnalysisInfo,"%s\n","_________________________________________________________");
  fprintf(AnalysisInfo,"%s\t\t %f\n","Minimum Distance:",min_distance);
  fprintf(AnalysisInfo,"%s\t\t %f\n","Maximum Distance:",max_distance);
  fprintf(AnalysisInfo,"%s\t %.2f %s\n","Fraction of Weighted Events:",fev_weighted,"%");
  fprintf(AnalysisInfo,"%s\n\n","---------------------------------------------------------");
  fclose(AnalysisInfo);
  
  cout<<"Minimum Distance: "<<min_distance<<endl;
  cout<<"Maximum Distance: "<<max_distance<<endl;
  //cout<<"Neighbors: "<<neighbors<<endl;
  
  //data->GetEntry(chosen_id);
  //for(int k=0; k<6; k++){
    //mc_chosen[k][0] = data_stored[k][0];
    //mc_chosen[k][1] = data_stored[k][1];
    //mc_chosen[k][2] = data_stored[k][2];
  //}
  //results->Fill();
  
  TFile *arq = new TFile("Results.root","recreate");
  results->Write();
  arq->Close();
  
}//End main program
