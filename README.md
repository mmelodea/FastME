________________________________________________________________________________________
<img align="left" src="https://github.com/mmelodea/FastMatrixElement/blob/master/FastMatrixElement/.fme_logo.png" width="180"> <h1>FastME</h1>
<h6>Probing and Discriminating Data by MC Topology (stable & updated version of the project)</h6>
________________________________________________________________________________________
FastME stands for Fast Matrix Element and its main goal is achieve similar capabability and full physic information as the many types of Matrix Element method available today. However, FastME try to reduce the computing time needed to get such info.  
The Monte Carlo events already contain the full physic information (that can be in higher orders of parameters correction) about the related event and FastME uses these events as inputs to compute in a fast way the nature (signal/background) of a specific event (a data event, for instance). The process relies into compare de available parameters (pT and eta) between the event and the MC. Such comparison allows one to compute the distance between an event and each MC in that phase space. Finally, the probability of an analised event be of one type or another type (specified by MC) is computed using the minimum distances found to each one of the MC used.


<h5>Usage and Setup</h5>  
This package is designed to be plugged into CMSSW. However, it also can be compiled by a c++ compiler if one has the ROOT6 libraries available and support to the new c++ version (c++11).  
Into CMSSW you just need to do 'source SetUpEnvironment.sh' (the script inside the package) and then all the needed configuration will be setup and the package will be compiled automatically. The executable generated is called "fastme" (present on CMSSW_BASE/bin). For checking the commands available use "fastme -help".  
To use fastme the user has to pass a configuration file that specifies some parameters in the analysis and also for control the results flow during/after the analysis. To get more info look in the file "fme_config.dat" to see the setup available in the package. Note that the root files input to the software need a specific format (as the example "fme_ntuple_format.root"). The package allows one to process multiple Data and MC files in just one run. Each Data event is compared to each MC event from all MC samples at the same time. The results from each Data file are identifiable by a tag based on its order in the config file. There are 2 methods to compare the Data and MC particles: minimum distance (mindr) and average of identical particles (mean).
The software produces essencially a ROOT TTree object containing usefull information used by itself. The tree named "fme_tree" is an internal class of info that is used in the further option in the software that makes the discriminant distribution analysis.

The 2 main options to use the software are:

<h6>fastme -a fme_config.dat</h6>
This option will run the software to find the closest MC events to each data event and it will store in the tree fme_tree the data file analised, the MC file matched to the data events, the minimum distances found.


<h6> fastme -d fme_config.dat</h6>
This option will take the previous file created and the tree stored in it and it will compute the discriminant acording to the scheme defined by the user (set of sig/bkg_data and sig/bkg_mc). It also will provide the discriminant variation for each MC file, which allows one to know which of the MC is more similar/different to the analised data.
