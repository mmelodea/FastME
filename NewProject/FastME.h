#ifndef FME_h
#define FME_h

#include <iostream>

class FME{
  
public:
  FME(){}
  ~FME(){}
  double compute_DR(string type, Float_t *Data, Float_t *MC);
  
private:
  double FS4l_DR(Float_t *Data, Float_t *MC);
  double FS4l2j_DR(Float_t *Data, Float_t *MC);
  double FSlv2j_DR(Float_t *Data, Float_t *MC);
};

#endif