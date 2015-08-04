///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::: 4 LEPTONS + 2 JETS FINAL STATE CASE ::::::::::::::::::::::::

#ifndef FS4l2j_DrComputers_h
#define FS4l2j_DrComputers_h

#include <TString.h>

double	FS4l2j_Res_DrOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution);
  double	FS4l2j_Res_DrOrder_Reso  (Float_t Data[6][3][2], Float_t MC[6][3][2]);
  double	FS4l2j_Res_DrOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]);

double	FS4l2j_noRes_DrOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution);
  double	FS4l2j_noRes_DrOrder_Reso  (Float_t Data[6][3][2], Float_t MC[6][3][2]);
  double	FS4l2j_noRes_DrOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]);

double	FS4l2j_Res_PtOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution);
  double	FS4l2j_Res_PtOrder_Reso  (Float_t Data[6][3][2], Float_t MC[6][3][2]);
  double	FS4l2j_Res_PtOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]);

double	FS4l2j_noRes_PtOrder(Float_t Data[6][3][2], Float_t MC[6][3][2], TString Resolution);
  double	FS4l2j_noRes_PtOrder_Reso  (Float_t Data[6][3][2], Float_t MC[6][3][2]);
  double	FS4l2j_noRes_PtOrder_noReso(Float_t Data[6][3][2], Float_t MC[6][3][2]);

#endif