#script to set up all the needed configuration and to compile FastME package

echo "Creating CMSSW_8_0_2 release..."
export SCRAM_ARCH=slc6_amd64_gcc493
scram project CMSSW CMSSW_8_0_2
cd CMSSW_8_0_2/src/

echo "Setting up CMSSW_8_0_2 environment..."
cmsenv

echo "Compiling FastME package..."
scram b -v 
