///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifndef FSlv2j_DrComputer_h
#define FSlv2j_DrComputer_h

#include <TString.h>

///FINAL STATE lvjj
double	FSlv2j_Jet_DrOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString Resolution);
  double	FSlv2j_jet_DrOrder_Reso(Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FSlv2j_jet_DrOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]);

double	FSlv2j_Jet_PtOrder(Float_t Data[4][3][2], Float_t MC[4][3][2], TString Resolution);
  double	FSlv2j_jet_PtOrder_Reso(Float_t Data[4][3][2], Float_t MC[4][3][2]);
  double	FSlv2j_jet_PtOrder_noReso(Float_t Data[4][3][2], Float_t MC[4][3][2]);

#endif