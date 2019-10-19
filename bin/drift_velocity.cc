
#include <iostream>
#include <math.h>

int main() {
  double efield = 0.5;
  double temperature = 89.2;


  double tshift = -87.203+temperature;
  double xFit = 0.0938163-0.0052563*tshift-0.0001470*tshift*tshift;
  double uFit = 5.18406+0.01448*tshift-0.003497*tshift*tshift-0.000516*tshift*tshift*tshift;
  double vd;

  double  P1 = -0.04640; // K^-1                                                                                                                                               
  double  P2 = 0.01712;  // K^-1                                                                                                                                               
  double  P3 = 1.88125;   // (kV/cm)^-1                                                                                                                                        
  double  P4 =  0.99408;    // kV/cm                                                                                                                                           
  double  P5 =  0.01172;   // (kV/cm)^-P6                                                                                                                                      
  double  P6 =  4.20214;
  double  T0 =  105.749;  // K                                                                                                                                                 

  double   P1W = -0.01481; // K^-1                                                                                                                                            
  double   P2W = -0.0075;  // K^-1                                                                                                                                              
  double   P3W =  0.141;   // (kV/cm)^-1                                                                                                                                       
  double   P4W =  12.4;    // kV/cm                                                                                                                                            
  double   P5W =  1.627;   // (kV/cm)^-P6                                                                                                                                      
  double   P6W =  0.317;
  double   T0W =  90.371;  // K                                                                                                                                                
  if (efield < xFit) vd=efield*uFit;
  else if (efield<0.619) {
    vd = ((P1*(temperature-T0)+1)
	  *(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
	  +P2*(temperature-T0));
  }
  else if (efield<0.699) {
    vd = 12.5*(efield-0.619)*((P1W*(temperature-T0W)+1)
			      *(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
			      +P2W*(temperature-T0W))+
      12.5*(0.699-efield)*((P1*(temperature-T0)+1)
			   *(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
			   +P2*(temperature-T0));
  }
  else {
    vd = ((P1W*(temperature-T0W)+1)
	  *(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
	  +P2W*(temperature-T0W));
  }
  
  vd /= 10.;
  
  std::cout<<vd<<std::endl;
  return 0;
}
