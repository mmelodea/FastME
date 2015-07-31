::::::::::::::::::: FastME ::::::::::::::::::::
 Probing and Discriminating Data Code MC Based
:::::::::::::::::::::::::::::::::::::::::::::::

Use Instructions:

1. You need format your ntuples to get a matrix with the ordered objects (electron, muon, jets, met ...) and properties (pT, eta and phi):
 
 1.1. The matrix must be tri-dimensional (matrix[number_final_state][objects][properties]):
  1.1.1. The first indice indicate the object order;
  1.1.2. The second indice indicate the properties order;
  1.1.3. The third indice indicate the propertie and its resolution order;
  
 Example- For the 4l (full leptonic case):
     Data[0][2][0] point to some lepton phi value;
     While, Data[0][2][1] point to lepton phi resolution value;
     
2. To run the analysis I recommend to use the interface code (at path "Interface/"). The code there make the full connection with
   the Fast Matrix Element library. The compilation was originaly processed on ROOT (I recommend to use it because the its libraries
   included in the FastME library). To use that code (there is a example file there):
   
   2.1. Write a txt file with data, sig, bkg ntuples path, the TTree name and the Branch name corresponding to matrix.
   This file must be passed as an argument to the function contained in the interfacing code;
   **(!!Do not forget to put "#fim" at the end of txt file!!)
   2.2. In ROOT command line make the argument function compile command, pointing to the code at "Interface/";
   2.3. After compilation, use the function instantiated to call FastME;
   
   Simplify Steps:  root [0].L runFME.C+
		    root [1]runFME("Inputs")
		    
3. The results of the analysis will be a ROOT file containing a Tree with quantities like minDR_ToSig/Bkg (minimum data-sig/bkg distances found),
   the P_SB discriminant value (indicating how much sig/bkg the data look like), the fraction (%) of event inputs like sig/bkg, and the event
   weight (there are 3 weight defined: one given by sig, one given by bkg and one given by the relation between the two first).