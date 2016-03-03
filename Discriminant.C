///############################################################################################################$
///======================[ Discriminant Module - Analysisng the Results from Distances Found ]==================
///=====================================[ Code Author: Miqueias M. de Almeida ]================================
///############################################################################################################$



#include <iostream>
#include <vector>
#include <TTree.h>
#include <TStopwatch.h>

TTree *Discriminant(TTree *mtree, UInt_t N_Cores, Int_t nData, const Int_t N_MCT, Int_t verbose){

  TStopwatch t2;
  t2.Start();
  ///_______________________ Compute discriminant from MDMCED _________________________________________________$
  cout<<":::::::                                                                                :::::::"<<endl;
  cout<<":::::::::::                                                                        :::::::::::"<<endl;
  cout<<"::::::::::::::::::                                                          ::::::::::::::::::"<<endl;
  cout<<"::::::::::::::::::::::::::::::::[ Computing discriminant ]::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<"::::::::::::::::::                                                          ::::::::::::::::::"<<endl;
  cout<<":::::::::::                                                                        :::::::::::"<<endl;
  cout<<":::::::                                                                                :::::::"<<endl;
  ///----------------------------------------------------------------------------------------------------------$

  ///Set the input tree
  Int_t iEvent, TMcType, Indice;
  Double_t Mdist;
  mtree->SetBranchAddress("iEvent",&iEvent);
  mtree->SetBranchAddress("Mdist",&Mdist);
  mtree->SetBranchAddress("TMcType",&TMcType);
  mtree->SetBranchAddress("Indice",&Indice);  
  Int_t fentries = mtree->GetEntries();

  vector<Int_t> McIndex, McCat;
  Double_t Global_PsbDist;
  vector<Double_t> MinDist, Local_PsbDist;
  TTree *ftree = new TTree("FastME","Fast Matrix Element Analysis Results");
  ftree->SetDirectory(0);
  ftree->Branch("McIndex","vector<Int_t>",&McIndex);
  ftree->Branch("McCat","vector<Int_t>",&McCat);
  ftree->Branch("MinDist","vector<Double_t>",&MinDist);
  ftree->Branch("Global_PsbDist",&Global_PsbDist,"Global_PsbDist/D");
  ftree->Branch("Local_PsbDist","vector<Double_t> Local_PsbDist",&Local_PsbDist);

  ///Find the tree sectors
  Int_t TreeSectors[N_Cores];
  if(fentries % N_Cores != 0){
    cout<<"[Error] Something gone wrong... non-integer tree sectors!!"<<endl;
    throw exception();
  }
  
  Int_t EndSector = fentries/N_Cores; //How many entries in each core job                                     
  for(Int_t ic=0; ic<(Int_t)N_Cores; ic++) TreeSectors[ic] = ic*EndSector;
  
  ///Getting results from analysis
  for(Int_t data=0; data<nData; data++){
    if( verbose != 0 && data%(nData/10) == 0){
      cout<< Form(":: [Remaining]:   %i Events",nData-data)<<"\t\t";
      t2.Print();
    }

    Int_t MinSigIndex = -99, GMinBkgIndex = -99;
    Double_t min_dr_sig = 1.e15, global_min_dr_bkg = 1.e15;
    Double_t local_min_dr_bkg[N_MCT], local_min_bkg_index[N_MCT];
    McIndex.clear();
    McCat.clear();
    MinDist.clear();
    Local_PsbDist.clear();
    for(Int_t p=0; p<N_MCT; p++){
      local_min_bkg_index[p] = -99;
      local_min_dr_bkg[p] = 1.e15;
      McCat.push_back( -99 );
      McIndex.push_back( -99 );
      MinDist.push_back( -99. );
      Local_PsbDist.push_back( -99. );
    }
    for(Int_t ic=0; ic<(Int_t)N_Cores; ic++){
      mtree->GetEntry(TreeSectors[ic]+data);//TreeSectors aligns the results from different cores
      if(iEvent != data){
	cout<<"[Error] Something gone wrong... iEvent != data"<<endl;
	throw exception();
      }
      
      ///Finds closet MC Signal
      if(TMcType == 0)
	if( Mdist < min_dr_sig ){
	  min_dr_sig = Mdist;
	  MinSigIndex = Indice;
	  McCat[0] = 0;
	}

      ///Finds closet MC Background
      if(TMcType > 0){
	///The general most close MC Background
        if( Mdist < global_min_dr_bkg ) global_min_dr_bkg = Mdist;
	///Each MC Background
	if( Mdist < local_min_dr_bkg[TMcType] ){
	  McCat[TMcType] = TMcType;
	  local_min_dr_bkg[TMcType] = Mdist;
	  local_min_bkg_index[TMcType] = Indice;
	}
      }
    }

    MinDist[0] = min_dr_sig;
    McIndex[0] = MinSigIndex;
    Global_PsbDist = PsbD(min_dr_sig, global_min_dr_bkg);

    for(Int_t im=1; im<N_MCT; im++){
      McIndex[im] = local_min_bkg_index[im];
      MinDist[im] = local_min_dr_bkg[im];
      Local_PsbDist[im] = PsbD(min_dr_sig, local_min_dr_bkg[im]);
    }
    if( verbose > 1 )
      cout<< Form("GSigMin:   %f\t\tGBkgMin:   %f\t\tGPsbDMinDist:   %f", min_dr_sig, global_min_dr_bkg, Global_PsbDist) << endl;
    
    ftree->Fill();
  }
  
  ///________________________________ Stoping timming ________________________________________________________
  cout<<"\n::::::::::::::::::::::::::::::::::::[ Process Finished ]::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: [Computing Total Time]: "; t2.Stop(); t2.Print();
  cout<<":: [Sending TTree results...]"<<endl;
  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"<<endl;
  ///---------------------------------------------------------------------------------------------------------

  return ftree;
}
