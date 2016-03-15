#script to set up all the needed configuration and to compile FastME package

echo "Creating CMSSW_8_0_2 release..."
cd ../
export SCRAM_ARCH=slc6_amd64_gcc493
scram project CMSSW CMSSW_8_0_2
mv FastMatrixElement CMSSW_8_0_2/src/
cd CMSSW_8_0_2/src/

echo "Setting up CMSSW_8_0_2 environment..."
cmsenv

echo "Compiling FastME package..."
scram b -v

echo "Proccess finished!"
