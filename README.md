________________________________________________________________________________________
<img align="left" src="https://raw.githubusercontent.com/mmelodea/FastME/master/fme_logo.png" width="180"> <h1>FastME</h1>
<h6>Probing and Discriminating Data by MC Topology (stable & updated version of the project)</h6>
________________________________________________________________________________________
FastME stands for Fast Matrix Element and its main goal is achieve similar capabability and full physic information as the many types of Matrix Element method available today. However, FastME try to reduce the computing time needed to get such info.  
The Monte Carlo events already contain the full physic information (that can be in higher orders of parameters correction) about the related event and FastME uses these events as inputs to compute in a fast way the nature (signal/background) of a specific event (a data event, for instance). The process relies into compare de available parameters (pT and eta) between the event and the MC. Such comparison allows one to compute the distance between an event and each MC in that phase space. Finally, the probability of an analised event be of one type or another type (specified by MC) is computed using the minimum distances found to each one of the MC used.


<h5>Usage and Setup</h5>  
This package is designed to be plugged into CMSSW. However, it also can be compiled by a c++ compiler since one has the ROOT6 libraries available and support to the new c++ version (c++11).  
Into CMSSW you just need to do 'source SetUpEnvironment.sh' and then all the needed configuration will be setup and the package will be compiled. The executable generated is called "fastme" (present on CMSSW_BASE/bin). For check the commands available use "fastme -help".


<h5>Results</h5>
The branches present in the resulting file store the follow informations:  
(1) Global_PsbDist - the discriminant computed taking in account all MCs (a minimum DR is chosen between the MC backgrounds);  
(2) Local_PsbDist - the discriminant computed to each of the MC backgrounds (its value is -99 for signal, since doesn't make sense compare signal to signal DR; also, note that the framework assumes only one signal source by analysis);  
(3) McCat - the MC category: (always!)0 for signal and >0 for backgrounds;  
(4) McIndex - the MC chosen as the most close to the analised event;  
(5) MinDist - the distance found for the most close MC found to the analised event.
