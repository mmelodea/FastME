///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS			::::::::::
///::::::		    Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::: LEPTONS + MET + 2 JETS FINAL STATE CASE :::::::::::::::::::::::::::

#include "FSlv2j_DrComputers.h"
#include <TString.h>

///JET DISTANCE ORDERING
double FSlv2j_Jet_DrOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString Resolution){
  int SetResolution = (Resolution == "True")? 0 : 1;
  switch(SetResolution){
    case 0:
      return FSlv2j_jet_DrOrder_Reso(Data,MC);
      break;
    case 1: 
      return FSlv2j_jet_DrOrder_noReso(Data,MC);
      break;
    default:
      return FSlv2j_jet_DrOrder_noReso(Data,MC);
      break;
  }
}

///USES RESOLUTION
double FSlv2j_jet_DrOrder_Reso(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start;

    ///Electron
    start = 0;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    ///MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    ///-------------------------------   Jets   --------------------------------------
    start = 2;
    double fdpt2    = pow( (Data[start][0][0]-MC[start][0][0])/Data[start][0][1] ,2 );
    double fdeta2   = pow( (Data[start][1][0]-MC[start][1][0])/Data[start][1][1] ,2 );
    double fdphi2   = pow( (Data[start][2][0]-MC[start][2][0])/Data[start][2][1] ,2 );
    double sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
    double  sdpt2    = pow( (Data[start][0][0]-MC[start+1][0][0])/Data[start][0][1] ,2 );
    double  sdeta2   = pow( (Data[start][1][0]-MC[start+1][1][0])/Data[start][1][1] ,2 );
    double  sdphi2   = pow( (Data[start][2][0]-MC[start+1][2][0])/Data[start][2][1] ,2 );
    double  sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
    if(sum_dr1 < sum_dr2){
      sum_dpt2  += fdpt2;
      sum_deta2 += fdeta2;
      sum_dphi2 += fdphi2;
  
      sum_dpt2  += pow( (Data[start+1][0][0]-MC[start+1][0][0])/Data[start+1][0][1] ,2 );
      sum_deta2 += pow( (Data[start+1][1][0]-MC[start+1][1][0])/Data[start+1][1][1] ,2 );
      sum_dphi2 += pow( (Data[start+1][2][0]-MC[start+1][2][0])/Data[start+1][2][1] ,2 );
    }
      
    if(sum_dr2 < sum_dr1){
      sum_dpt2  += sdpt2;
      sum_deta2 += sdeta2;
      sum_dphi2 += sdphi2;

      sum_dpt2  += pow( (Data[start+1][0][0]-MC[start][0][0])/Data[start+1][0][1] ,2 );
      sum_deta2 += pow( (Data[start+1][1][0]-MC[start][1][0])/Data[start+1][1][1] ,2 );
      sum_dphi2 += pow( (Data[start+1][2][0]-MC[start][2][0])/Data[start+1][2][1] ,2 );
    }

  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}

///NO USES RESOLUTION
double FSlv2j_jet_DrOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start;

    ///Electron
    start = 0;
    dpt  = (Data[start][0][0] - MC[start][0][0]);
    deta = (Data[start][1][0] - MC[start][1][0]);
    dphi = (Data[start][2][0] - MC[start][2][0]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    ///MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0]);
    deta = (Data[start][1][0] - MC[start][1][0]);
    dphi = (Data[start][2][0] - MC[start][2][0]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    ///-------------------------------   Jets   --------------------------------------
    start = 2;
    double fdpt2    = pow( Data[start][0][0]-MC[start][0][0] ,2 );
    double fdeta2   = pow( Data[start][1][0]-MC[start][1][0] ,2 );
    double fdphi2   = pow( Data[start][2][0]-MC[start][2][0] ,2 );
    double sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
    double  sdpt2    = pow( Data[start][0][0]-MC[start+1][0][0] ,2 );
    double  sdeta2   = pow( Data[start][1][0]-MC[start+1][1][0] ,2 );
    double  sdphi2   = pow( Data[start][2][0]-MC[start+1][2][0] ,2 );
    double  sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
    if(sum_dr1 < sum_dr2){
      sum_dpt2  += fdpt2;
      sum_deta2 += fdeta2;
      sum_dphi2 += fdphi2;
  
      sum_dpt2  += pow( Data[start+1][0][0]-MC[start+1][0][0] ,2 );
      sum_deta2 += pow( Data[start+1][1][0]-MC[start+1][1][0] ,2 );
      sum_dphi2 += pow( Data[start+1][2][0]-MC[start+1][2][0] ,2 );
    }
      
    if(sum_dr2 < sum_dr1){
      sum_dpt2  += sdpt2;
      sum_deta2 += sdeta2;
      sum_dphi2 += sdphi2;

      sum_dpt2  += pow( Data[start+1][0][0]-MC[start][0][0] ,2 );
      sum_deta2 += pow( Data[start+1][1][0]-MC[start][1][0] ,2 );
      sum_dphi2 += pow( Data[start+1][2][0]-MC[start][2][0] ,2 );
    }

  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///==============================================================================================



///JETS PT ORDERING
double FSlv2j_Jet_PtOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString Resolution){
  int SetResolution = (Resolution == "True")? 0 : 1;
  switch(SetResolution){
    case 0:
      return FSlv2j_jet_PtOrder_Reso(Data,MC);
      break;
    case 1: 
      return FSlv2j_jet_PtOrder_noReso(Data,MC);
      break;
    default:
      return FSlv2j_jet_PtOrder_noReso(Data,MC);
      break;
  }
}

///USES RESOLUTION
double FSlv2j_jet_PtOrder_Reso(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start;

    ///Electron
    start = 0;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    ///MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0])/Data[start][0][1];
    deta = (Data[start][1][0] - MC[start][1][0])/Data[start][1][1];
    dphi = (Data[start][2][0] - MC[start][2][0])/Data[start][2][1];
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    ///-----------------------------   Jets   --------------------------------------
    start = 2;
    int data_st = (Data[start][0][0] > Data[start+1][0][0])? start : start+1;
    int data_nd = (Data[start][0][0] < Data[start+1][0][0])? start : start+1;
    int   mc_st = (  MC[start][0][0] > MC[start+1][0][0]  )? start : start+1;
    int   mc_nd = (  MC[start][0][0] < MC[start+1][0][0]  )? start : start+1;
          
    //The Highest pT jets
    sum_dpt2  += pow( (Data[data_st][0][0]-MC[mc_st][0][0])/Data[data_st][0][1] ,2 );
    sum_deta2 += pow( (Data[data_st][1][0]-MC[mc_st][1][0])/Data[data_st][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_st][2][0]-MC[mc_st][2][0])/Data[data_st][2][1] ,2 );

    //The last jets (smaller pT)
    sum_dpt2  += pow( (Data[data_nd][0][0]-MC[mc_nd][0][0])/Data[data_nd][0][1] ,2 );
    sum_deta2 += pow( (Data[data_nd][1][0]-MC[mc_nd][1][0])/Data[data_nd][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_nd][2][0]-MC[mc_nd][2][0])/Data[data_nd][2][1] ,2 );
    
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}

///USES RESOLUTION
double FSlv2j_jet_PtOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]){
  double dpt = 0, deta = 0, dphi = 0;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start;

    ///Electron
    start = 0;
    dpt  = (Data[start][0][0] - MC[start][0][0]);
    deta = (Data[start][1][0] - MC[start][1][0]);
    dphi = (Data[start][2][0] - MC[start][2][0]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
    
    ///MET
    start = 1;
    dpt  = (Data[start][0][0] - MC[start][0][0]);
    deta = (Data[start][1][0] - MC[start][1][0]);
    dphi = (Data[start][2][0] - MC[start][2][0]);
    sum_dpt2  += pow(dpt,2);
    sum_deta2 += pow(deta,2);
    sum_dphi2 += pow(dphi,2);
  
    ///-----------------------------   Jets   --------------------------------------
    start = 2;
    int data_st = (Data[start][0][0] > Data[start+1][0][0])? start : start+1;
    int data_nd = (Data[start][0][0] < Data[start+1][0][0])? start : start+1;
    int   mc_st = (  MC[start][0][0] > MC[start+1][0][0]  )? start : start+1;
    int   mc_nd = (  MC[start][0][0] < MC[start+1][0][0]  )? start : start+1;
          
    //The Highest pT jets
    sum_dpt2  += pow( Data[data_st][0][0]-MC[mc_st][0][0] ,2 );
    sum_deta2 += pow( Data[data_st][1][0]-MC[mc_st][1][0] ,2 );
    sum_dphi2 += pow( Data[data_st][2][0]-MC[mc_st][2][0] ,2 );

    //The last jets (smaller pT)
    sum_dpt2  += pow( Data[data_nd][0][0]-MC[mc_nd][0][0] ,2 );
    sum_deta2 += pow( Data[data_nd][1][0]-MC[mc_nd][1][0] ,2 );
    sum_dphi2 += pow( Data[data_nd][2][0]-MC[mc_nd][2][0] ,2 );
    
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///====================================================================================================