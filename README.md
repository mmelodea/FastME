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
It's still possible filter Data events with more particle than the expected for the final state. Remembering that always the MC samples are templates that must contain only the final particles. Be aware that for 'mean' method is not possible to do that, since the software doesn't find the closest particles.


<h5>Results</h5>
The software produces 2 ROOT TTrees objects containing usefull information used by itself. One of them "fme_tree" is an internal class of info that by default is not stored in the final ROOT file, but it can be usefull depending on some study. The TTree "FastME" is the final ROOT file that contains all the important results for the user and has the following branches:

(1) Global_PsbDist - the discriminant computed taking into account all MCs (a minimum DR is chosen between the MC backgrounds);  
(2) Local_PsbDist - the discriminant computed to each of the MC backgrounds (its value is -99 for signal, since doesn't make sense compare signal to signal DR; also, note that the framework assumes only one signal source by analysis);
(3) PairedMC - the closest MC event found, for each MC sample, to the analised event;
(4) PairedMCType - the MC category: (always!)0 for signal and >0 for backgrounds;  (be carefull when writing the MC files in the config file)
(5) MinDistance - the distance computed for the closest MC event found.


<h5>View from the Software Running</h5>
<img align="center" src="https://github.com/mmelodea/FastMatrixElement/blob/master/FastMatrixElement/.SoftwareRunning.png" width="800">


<h5>Plots that can be done through the software for further studies</h5>
The plots below show the main results from the analysis. From left to right, up to bottom, the plots are: the discriminant distribution (and the respective signal significance) for all signal and background test files, the minimum distances to the paired MCs, the indexes of the paired MCs to signal test evens, the ROC curve based on the discriminants, the discriminant distribution for each data test and again the indexes of paired MCs (now to background test events)
<img align="center" src="https://github.com/mmelodea/FastMatrixElement/blob/master/FastMatrixElement/.fme_new_format.png" width="1000">
