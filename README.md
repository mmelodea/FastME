# FastME
<img src="https://raw.githubusercontent.com/mmelodea/FastME/master/fme_logo.png" width="200">

Probing and Discriminating Data by MC Topology
  (stable & updated version of the project)

FastME stands for Fast Matrix Element and its main goal is achieve similar capabability and full physic information as the many types of Matrix Element method available today. However, FastME try to reduce the computing time needed to get such info.  
The Monte Carlo events already contain the full physic information about the related event and FastME use these events as inputs to compute in a fast way the nature of a specific event (a data event, for instance). The process relies into compare de available parameters (pT and eta) between the event and the MC. Such comparison allows one to compute the distance between an event and each MC in that phase space. Finally, the probability of an analised event be of one type or another type (specified by MC) is computed using those distances.