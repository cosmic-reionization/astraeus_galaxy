#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

//#include "const.h"
//#include "utils.h"

//#include "dconfObj.h"
//#include "gal_gtree.h"
#include "starformation_SNfeedback.h"

/*---------------------------------------------------------------------*/
/* SUPERNOVA FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

float calc_SNejection_efficiency(const int delayedSNfeedback, float &FW, float *SNenergy, const int32_t &snapnumber, float &Mvir, float &Rvir, const float &scalefactor, float &fracMgasMer, float &MgasIni, float corrFactorTimeStep){
  // Graviational constant kpc Msol^-1 (km/s)^2
  float G = 4.302e-6; 

  // Free-fall velocity??? in km/s
  float vc = sqrt(G * Mvir / (Rvir * scalefactor));

  //   printf("Mvir = %e\t rvir = %e\t vc = %e\n", thisGal->Mvir, thisGal->Rvir, vc);
  float fw = FW / (1. + (corrFactorTimeStep - 1.) * fracMgasMer);
  float SNejection_efficiency = 0.;
  float Mgas = MgasIni;

  if(delayedSNfeedback == 1){
    float SNenergyPast = calc_SNenergy_past(snapnumber, SNenergy);
    float factorSNenergyPast = 0.;

    if(Mgas > 0.){
      factorSNenergyPast = 1. - fw * SNenergyPast / (Mgas * vc * vc);
    }
    
    float SNenergyCurrent = get_SNenergy_current(snapnumber, SNenergy);
    float factorSNenergyCurrent = vc * vc / (vc * vc + fw * SNenergyCurrent);
    
    SNejection_efficiency = factorSNenergyCurrent * fmax(0., factorSNenergyPast);
  }
  else {
    // Fixed free-fall velocity in km/s
    float vs = 611.;
    SNejection_efficiency = vc * vc / (vc * vc + fw * vs * vs);
  }
  
  return SNejection_efficiency;
}


float calc_remaining_gasfraction(const int delayedSNfeedback, float &FW, float *SNenergy, const int32_t &snapnumber, float &Mvir, float &Rvir, const float &scalefactor, float &fracMgasMer, float &MgasIni, float corrFactorTimeStep, float &fej, float &feff){
  // Graviational constant kpc Msol^-1 (km/s)^2
  float G = 4.302e-6; 

  // Circular velocity in km/s
  float vc = sqrt(G * Mvir / (Rvir * scalefactor));

  float fw = FW / (1. + (corrFactorTimeStep - 1.) * fracMgasMer);
  float remainingGasfraction = 0.;
  
  if(delayedSNfeedback == 1)
  {
    float SNenergyCurrent = get_SNenergy_current(snapnumber, SNenergy);
    
    remainingGasfraction = (fej - feff) * (vc * vc + fw * SNenergyCurrent) / (vc * vc * (1. - feff));
  }
  else{
    remainingGasfraction = feff * (1. - fej) / (fej * (1. - feff));
  }
  
  return remainingGasfraction;
}

float get_SNenergy_current(const int &snapnumber, float *SNenergy){
  float result = 0.;
  int snap = snapnumber;

  if(snap > 0){
    result = SNenergy[snap*(snap-1)/2 + snap - 1];
  }

  return result;
}

float calc_SNenergy_past(const int &snapnumber, float *SNenergy){
  float result = 0.;
  int snap = snapnumber;
  //TODO: we provably need to change this. Currently, it seems it's accessing the memory of the SFR of the same halos in the previous snapshots, but in this implementation, we cannot do this. Consider passing the comulative stellar mass formed in the previous snapshots as an argument to this function and evolve_gal.
  //TODO: see Equation (8) in paper I https://arxiv.org/pdf/2004.08401
  for(int prevSnap=1; prevSnap<snap; prevSnap++){
    result += SNenergy[snap*(snap-1)/2 + prevSnap - 1]; //TODO:* thisGal->stellarmasshistory[prevSnap];
  }
  
  return result;
}


