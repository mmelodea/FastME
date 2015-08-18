///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::: 4 LEPTONS + 2 JETS FINAL STATE CASE ::::::::::::::::::::::::

#include "FS4l2j_DrComputers.h"
#include <TString.h>

///DISTANCE ORDERING + RESSONANCE CRITERY
double FS4l2j_Res_DrOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution){
  int SetResolution = (Resolution == "True")? 0 : 1;
  switch(SetResolution){
    case 0:
      return FS4l2j_Res_DrOrder_Reso(Data,MC);
      break;
    case 1: 
      return FS4l2j_Res_DrOrder_noReso(Data,MC);
      break;
    default:
      return FS4l2j_Res_DrOrder_Reso(Data,MC);
      break;
  }
}

///USES RESOLUTION
double FS4l2j_Res_DrOrder_Reso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  do{
      sum_dr1 = 0;
      sum_dr2 = 0;
      
      fdpt2    = pow( (Data[start][0][0]-MC[start][0][0])/Data[start][0][1] ,2 );
      fdeta2   = pow( (Data[start][1][0]-MC[start][1][0])/Data[start][1][1] ,2 );
      fdphi2   = pow( (Data[start][2][0]-MC[start][2][0])/Data[start][2][1] ,2 );
      sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
      sdpt2    = pow( (Data[start][0][0]-MC[start+1][0][0])/Data[start][0][1] ,2 );
      sdeta2   = pow( (Data[start][1][0]-MC[start+1][1][0])/Data[start][1][1] ,2 );
      sdphi2   = pow( (Data[start][2][0]-MC[start+1][2][0])/Data[start][2][1] ,2 );
      sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
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

    start += 2; 
  }while(start < 6);
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  
  return event_distance;
}

///NO USES RESOLUTION
double FS4l2j_Res_DrOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double fdpt2, fdeta2, fdphi2, sdpt2, sdeta2, sdphi2, sum_dr1, sum_dr2;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int start = 0;
  
  do{
      sum_dr1 = 0;
      sum_dr2 = 0;
      
      fdpt2    = pow( Data[start][0][0]-MC[start][0][0] ,2 );
      fdeta2   = pow( Data[start][1][0]-MC[start][1][0] ,2 );
      fdphi2   = pow( Data[start][2][0]-MC[start][2][0] ,2 );
      sum_dr1 = sqrt(fdpt2 + fdeta2 + fdphi2);
      
      sdpt2    = pow( Data[start][0][0]-MC[start+1][0][0] ,2 );
      sdeta2   = pow( Data[start][1][0]-MC[start+1][1][0] ,2 );
      sdphi2   = pow( Data[start][2][0]-MC[start+1][2][0] ,2 );
      sum_dr2 = sqrt(sdpt2 + sdeta2 + sdphi2);
      
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

    start += 2; 
  }while(start < 6);
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  
  return event_distance;
}
///=============================================================================================


///DISTANCE ORDERING + NO RESSONANCE CRITERY
double FS4l2j_noRes_DrOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution){
  int SetResolution = (Resolution == "True")? 0 : 1;
  switch(SetResolution){
    case 0:
      return FS4l2j_noRes_DrOrder_Reso(Data,MC);
      break;
    case 1: 
      return FS4l2j_noRes_DrOrder_noReso(Data,MC);
      break;
    default:
      return FS4l2j_noRes_DrOrder_noReso(Data,MC);
      break;
  }
}

///USES RESOLUTION
double FS4l2j_noRes_DrOrder_Reso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double dpt2, deta2, dphi2, sum_dr = 0, limiar;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int chosen = -1, p[] = {-1,-1,-1,-1,-1,-1};
  double vchosen[] = {-1,-1,-1};
  
  //Choose the l/j1-Data most close to l/j1-MC
  for(int i=0; i<6; i++){
    limiar = 1.E15;
    for(int j=0; j<6; j++){
      if((j > 3 && i < 4) || (j < 4 && i > 3)) continue;       		                ///Avoid leptons/jets mixture
      if(j == p[0] || j == p[1] || j == p[2] || j == p[3] || j == p[4]) continue;	///Avoid recounting
      
      dpt2  = pow( (Data[i][0][0]-MC[j][0][0])/Data[i][0][1] ,2 );
      deta2 = pow( (Data[i][1][0]-MC[j][1][0])/Data[i][1][1] ,2 );
      dphi2 = pow( (Data[i][2][0]-MC[j][2][0])/Data[i][2][1] ,2 );
      sum_dr = sqrt(dpt2 + deta2 + dphi2);
      
      if(sum_dr < limiar){
	limiar     = sum_dr;
	chosen     = j;
	vchosen[0] = dpt2;
	vchosen[1] = deta2;
	vchosen[2] = dphi2;
      }
    }
    p[i] = chosen;
    
    //Makes the sum of deltapT2, deta2 and dphi2 to get the events distance
    sum_dpt2  += vchosen[0];
    sum_deta2 += vchosen[1];
    sum_dphi2 += vchosen[2];
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}

///NO USES RESOLUTION
double FS4l2j_noRes_DrOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double dpt2, deta2, dphi2, sum_dr = 0, limiar;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int chosen = -1, p[] = {-1,-1,-1,-1,-1,-1};
  double vchosen[] = {-1,-1,-1};
  
  //Choose the l/j1-Data most close to l/j1-MC
  for(int i=0; i<6; i++){
    limiar = 1.E15;
    for(int j=0; j<6; j++){
      if((j > 3 && i < 4) || (j < 4 && i > 3)) continue;       		                ///Avoid leptons/jets mixture
      if(j == p[0] || j == p[1] || j == p[2] || j == p[3] || j == p[4]) continue;	///Avoid recounting
      
      dpt2  = pow( Data[i][0][0]-MC[j][0][0] ,2 );
      deta2 = pow( Data[i][1][0]-MC[j][1][0] ,2 );
      dphi2 = pow( Data[i][2][0]-MC[j][2][0] ,2 );
      sum_dr = sqrt(dpt2 + deta2 + dphi2);
      
      if(sum_dr < limiar){
	limiar     = sum_dr;
	chosen     = j;
	vchosen[0] = dpt2;
	vchosen[1] = deta2;
	vchosen[2] = dphi2;
      }
    }
    p[i] = chosen;
    
    //Makes the sum of deltapT2, deta2 and dphi2 to get the events distance
    sum_dpt2  += vchosen[0];
    sum_deta2 += vchosen[1];
    sum_dphi2 += vchosen[2];
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///========================================================================================



///PT ORDERING + RESSONANCE CRITERY
double FS4l2j_Res_PtOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution){
  int SetResolution = (Resolution == "True")? 0 : 1;
  switch(SetResolution){
    case 0:
      return FS4l2j_Res_PtOrder_Reso(Data,MC);
      break;
    case 1: 
      return FS4l2j_Res_PtOrder_noReso(Data,MC);
      break;
    default:
      return FS4l2j_Res_PtOrder_noReso(Data,MC);
      break;
  }
}

///USES RESOLUTION
double FS4l2j_Res_PtOrder_Reso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double data_limiar, mc_limiar;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int data_ch = -1, mc_ch = -1, data_p[] = {-1,-1,-1,-1,-1,-1}, mc_p[] = {-1,-1,-1,-1,-1,-1};
  
  //Organizing leptons by pT
  for(int i=0; i<6; i++){
    data_limiar = 0.;
    mc_limiar   = 0.;
        
    for(int j=0; j<6; j++){
      if((j > 3 && i < 4) || (j < 4 && i > 3)) continue;       		        ///Avoid leptons/jets mixture
      
      if(j != data_p[0] && j != data_p[1] && j != data_p[2] && j != data_p[4]){         ///Avoid Data recounting
	if(Data[j][0][0] > data_limiar){
	  data_limiar = Data[j][0][0];
	  data_ch = j;
	}
      }
      if(j != mc_p[0] && j != mc_p[1] && j != mc_p[2] && j != mc_p[4]){			///Avoid MC recounting
	if(MC[j][0][0] > mc_limiar){
	  mc_limiar = MC[j][0][0];
	  mc_ch = j;
	}
      }
    }
    data_p[i] = data_ch;
    mc_p[i] = mc_ch;

    //Makes the sum over all final state particles
    sum_dpt2  += pow( (Data[data_ch][0][0]-MC[mc_ch][0][0])/Data[data_ch][0][1] ,2 );
    sum_deta2 += pow( (Data[data_ch][1][0]-MC[mc_ch][1][0])/Data[data_ch][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_ch][2][0]-MC[mc_ch][2][0])/Data[data_ch][2][1] ,2 );
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}

///NO USES RESOLUTION
double FS4l2j_Res_PtOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double data_limiar, mc_limiar;
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int data_ch = -1, mc_ch = -1, data_p[] = {-1,-1,-1,-1,-1,-1}, mc_p[] = {-1,-1,-1,-1,-1,-1};
  
  //Organizing leptons by pT
  for(int i=0; i<6; i++){
    data_limiar = 0.;
    mc_limiar   = 0.;
        
    for(int j=0; j<6; j++){
      if((j > 3 && i < 4) || (j < 4 && i > 3)) continue;       		        ///Avoid leptons/jets mixture
      
      if(j != data_p[0] && j != data_p[1] && j != data_p[2] && j != data_p[4]){         ///Avoid Data recounting
	if(Data[j][0][0] > data_limiar){
	  data_limiar = Data[j][0][0];
	  data_ch = j;
	}
      }
      if(j != mc_p[0] && j != mc_p[1] && j != mc_p[2] && j != mc_p[4]){			///Avoid MC recounting
	if(MC[j][0][0] > mc_limiar){
	  mc_limiar = MC[j][0][0];
	  mc_ch = j;
	}
      }
    }
    data_p[i] = data_ch;
    mc_p[i] = mc_ch;

    //Makes the sum over all final state particles
    sum_dpt2  += pow( Data[data_ch][0][0]-MC[mc_ch][0][0] ,2 );
    sum_deta2 += pow( Data[data_ch][1][0]-MC[mc_ch][1][0] ,2 );
    sum_dphi2 += pow( Data[data_ch][2][0]-MC[mc_ch][2][0] ,2 );
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///=====================================================================================



///PT ORDERING + NO RESSONANCE CRITERY
double FS4l2j_noRes_PtOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution){
  int SetResolution = (Resolution == "True")? 0 : 1;
  switch(SetResolution){
    case 0:
      return FS4l2j_noRes_PtOrder_Reso(Data,MC);
      break;
    case 1: 
      return FS4l2j_noRes_PtOrder_noReso(Data,MC);
      break;
    default:
      return FS4l2j_noRes_PtOrder_noReso(Data,MC);
      break;
  }
}

///USES RESOLUTION
double FS4l2j_noRes_PtOrder_Reso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int data_st = -1, data_nd = -1, mc_st = -1, mc_nd = -1, start = 0;
 
  //Organizing leptons by pT
  for(int i=0; i<3; i++){
    
    //Organizing data
    if(Data[start][0][0] > Data[start+1][0][0]){ data_st = start; data_nd = start+1; }
    if(Data[start][0][0] < Data[start+1][0][0]){ data_st = start+1; data_nd = start; }
      
    //Organizing MC
    if(MC[start][0][0] > MC[start+1][0][0]){ mc_st = start; mc_nd = start+1; }
    if(MC[start][0][0] < MC[start+1][0][0]){ mc_st = start+1; mc_nd = start; }

    //Makes the sum over all final state particles
    sum_dpt2  += pow( (Data[data_st][0][0]-MC[mc_st][0][0])/Data[data_st][0][1] ,2 );
    sum_deta2 += pow( (Data[data_st][1][0]-MC[mc_st][1][0])/Data[data_st][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_st][2][0]-MC[mc_st][2][0])/Data[data_st][2][1] ,2 );
    
    sum_dpt2  += pow( (Data[data_nd][0][0]-MC[mc_nd][0][0])/Data[data_nd][0][1] ,2 );
    sum_deta2 += pow( (Data[data_nd][1][0]-MC[mc_nd][1][0])/Data[data_nd][1][1] ,2 );
    sum_dphi2 += pow( (Data[data_nd][2][0]-MC[mc_nd][2][0])/Data[data_nd][2][1] ,2 );
    
    start += 2;
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}

///NO USES RESOLUTION
double FS4l2j_noRes_PtOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]){
  double sum_dpt2 = 0, sum_deta2 = 0, sum_dphi2 = 0, event_distance = -1;
  int data_st = -1, data_nd = -1, mc_st = -1, mc_nd = -1, start = 0;
 
  //Organizing leptons by pT
  for(int i=0; i<3; i++){
    
    //Organizing data
    if(Data[start][0][0] > Data[start+1][0][0]){ data_st = start; data_nd = start+1; }
    if(Data[start][0][0] < Data[start+1][0][0]){ data_st = start+1; data_nd = start; }
      
    //Organizing MC
    if(MC[start][0][0] > MC[start+1][0][0]){ mc_st = start; mc_nd = start+1; }
    if(MC[start][0][0] < MC[start+1][0][0]){ mc_st = start+1; mc_nd = start; }

    //Makes the sum over all final state particles
    sum_dpt2  += pow( Data[data_st][0][0]-MC[mc_st][0][0] ,2 );
    sum_deta2 += pow( Data[data_st][1][0]-MC[mc_st][1][0] ,2 );
    sum_dphi2 += pow( Data[data_st][2][0]-MC[mc_st][2][0] ,2 );
    
    sum_dpt2  += pow( Data[data_nd][0][0]-MC[mc_nd][0][0] ,2 );
    sum_deta2 += pow( Data[data_nd][1][0]-MC[mc_nd][1][0] ,2 );
    sum_dphi2 += pow( Data[data_nd][2][0]-MC[mc_nd][2][0] ,2 );
    
    start += 2;
  }
  
  event_distance = sqrt(sum_dpt2 + sum_deta2 + sum_dphi2);
  return event_distance;
}
///==========================================================================================
