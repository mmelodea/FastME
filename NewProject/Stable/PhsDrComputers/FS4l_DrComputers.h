///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifndef FS4l_DrComputers_h
#define FS4l_DrComputers_h

#include <TString.h>

///FINAL STATE 4 LEPTONS
double	FS4l_Res_DrOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString resolution);
  double	FS4l_Res_DrOrder_Reso  (Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FS4l_Res_DrOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]);
  
double	FS4l_noRes_DrOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString resolution);
  double	FS4l_noRes_DrOrder_Reso  (Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FS4l_noRes_DrOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]);

double	FS4l_Res_PtOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString resolution);
  double	FS4l_Res_PtOrder_Reso  (Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FS4l_Res_PtOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]);

double	FS4l_noRes_PtOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString resolution);
  double	FS4l_noRes_PtOrder_Reso  (Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FS4l_noRes_PtOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]);

#endif