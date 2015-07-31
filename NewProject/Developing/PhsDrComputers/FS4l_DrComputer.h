///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifndef FS4l_DrComputer_h
#define FS4l_DrComputer_h

///FINAL STATE 4 LEPTONS
double	FS4l_Res_DrOrder  (Float_t Data[4][3][2], Float_t MC[4][3][2]);
double	FS4l_noRes_DrOrder(Float_t Data[4][3][2], Float_t MC[4][3][2]);
double	FS4l_Res_PtOrder  (Float_t Data[4][3][2], Float_t MC[4][3][2]);
double	FS4l_noRes_PtOrder(Float_t Data[4][3][2], Float_t MC[4][3][2]);

#endif