///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::		PHASE ESPACE DATA-MC DISTANCE COMPUTERS	LIBRARY		::::::::::
///::::::		       Author: Miquéias M. de Almeida			::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifndef FS4l2j_DrComputer_h
#define FS4l2j_DrComputer_h

///FINAL STATE 4 LEPTONS
double	FS4l2j_Res_DrOrder  (Float_t Data[6][3][2], Float_t MC[6][3][2]);
double	FS4l2j_noRes_DrOrder(Float_t Data[6][3][2], Float_t MC[6][3][2]);
double	FS4l2j_Res_PtOrder  (Float_t Data[6][3][2], Float_t MC[6][3][2]);
double	FS4l2j_noRes_PtOrder(Float_t Data[6][3][2], Float_t MC[6][3][2]);

#endif