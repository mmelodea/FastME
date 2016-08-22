#ifndef InitScreen_h
#define InitScreen_h


#include <iostream>                                                            

void InitScreen(void){

  ///Fast Matrix Element screen openning
  std::cout<<std::endl;
  std::cout<<"                                                "<<ansi_red<<"   **             **     **************"<<ansi_reset<<std::endl;
  std::cout<<"                                                "<<ansi_red<<"  ***             ***   **************"<<ansi_reset<<std::endl;
  std::cout<<"       *********  ****     *********************"<<ansi_red<<" *****           *****  ***"<<ansi_reset<<std::endl;
  std::cout<<"      ********  ******   *********************  "<<ansi_red<<"*** ***         *** ***  ***"<<ansi_reset<<std::endl;
  std::cout<<"     ***      ***  ***  ***          ***        "<<ansi_red<<"***  ***       ***  ***   ************"<<ansi_reset<<std::endl;
  std::cout<<"    *****************************   ***         "<<ansi_red<<"***   ***     ***   ***   ************"<<ansi_reset<<std::endl;
  std::cout<<"   ******************************  ***          "<<ansi_red<<"***    ***   ***    ***  ***"<<ansi_reset<<std::endl;
  std::cout<<"  ***     ***    ***         ***  ***           "<<ansi_red<<"***     *** ***     *** ***"<<ansi_reset<<std::endl;
  std::cout<<" ***     ***    ***   *********  ***            "<<ansi_red<<" **       ***       **  **************"<<ansi_reset<<std::endl;
  std::cout<<"***     ***    ***  *********   ***             "<<ansi_red<<"  *                 *    **************"<<ansi_reset<<std::endl;
  std::cout<<std::endl;
  
  
  return;
}


#endif