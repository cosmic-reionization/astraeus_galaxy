#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#include "radiative_feedback.h"

/*---------------------------------------------------------------------*/
/* RADIATIVE FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

float calc_radiative_feedback(const float hubble_h, const float omega_m, const int32_t snapnumber, float Mvir, const float photHI_bg, const float scalefactor, const float zreion, const int32_t reion_flag, const float temp, const float mu, const int32_t radfeedback_model){
    float z = 1./scalefactor - 1.;
    //float zreion = thisGal->zreion;
    float fg = 1.;

    // TODO: remove the reion_flag and keep Gnedin as default
    if(reion_flag == 1){   
        if(zreion > 0.){
            if(radfeedback_model == 1){
                // requires input in units of Msun
                fg = fg_radfb_Sobacchi(photHI_bg, Mvir/hubble_h, z, zreion);
            }else if(radfeedback_model == 2){
                // requires inputs in units of h^-1 Msun
                fg = fg_radfb_Gnedin(Mvir, temp, mu, z, zreion, zreion, omega_m);
            }else if(radfeedback_model == 3){
                // requires inputs in units of h^-1 Msun
                fg = fg_radfb_Mjeans(Mvir, temp, mu, z, omega_m);
            }else if(radfeedback_model == 4){
                // requires inputs in units of h^-1 Msun
                fg = fg_radfb_Mcool(Mvir, temp, mu, z, omega_m);
            }else if(radfeedback_model == 5){
                // requires inputs in units of h^-1 Msun
                fg = fg_radfb_TempEvol(Mvir, temp, mu, z, zreion, zreion, omega_m);
            }
        } //else: Gnedin below
    }
    else{
        // requires inputs in units of h^-1 Msun
        fg = fg_radfb_Gnedin(Mvir, temp, mu, z, 8., 7., omega_m);
    }

    return fg;
}


/*---------------------------------------------------------------------*/
/* MASS LIMITS (all given in h^-1 Msun)                                */
/*---------------------------------------------------------------------*/

float get_Mcool(float Temp, float mu, float z, float omega_m){
    float sqrt_fact = sqrt(Temp*1.e-4/ (mu*(1.+z)));
    float Mcool = 5.24e8 * sqrt_fact * sqrt_fact * sqrt_fact / sqrt(omega_m);

    return Mcool;
}

float get_Mjeans(float Temp, float mu, float z, float omega_m){
    float sqrt_fact = sqrt(Temp*1.e-4/ (mu*(1.+z)));
    float Mjeans = 0.125 * 2.5e11 * sqrt_fact * sqrt_fact * sqrt_fact / sqrt(omega_m);

    return Mjeans;
}

float get_Mjeans_vir(float Temp, float mu, float z, float omega_m){
    float Mjeans = 0.075*get_Mjeans(Temp, mu, z, omega_m);

    return Mjeans;
}

/*---------------------------------------------------------------------*/
/* EXTREME RADIATIVE FEEDBACK DESCRIPTION                              */
/*---------------------------------------------------------------------*/
//TODO: it would be nice to give a reference to the paper that describes this model
float fg_radfb_Mjeans(float Mvir, float Temp, float mu, float z, float omega_m){
    float fg = 1.;

    float Mjeans = get_Mjeans_vir(Temp, mu, z, omega_m);

    fg = pow(2., -Mjeans/Mvir);

    return fg;
}

float fg_radfb_Mcool(float Mvir, float Temp, float mu, float z, float omega_m){
    float fg = 1.;

    float Mcool = get_Mcool(Temp, mu, z, omega_m);

    if(Mvir <= Mcool){
        fg = 0.;
    }

    return fg;
}


/*---------------------------------------------------------------------*/
/* RADIATIVE FEEDBACK DESCRIPTION FOLLOWING SOBACCHI 2013              */
/*---------------------------------------------------------------------*/
float fg_radfb_Sobacchi(float photHI, float Mvir, float z, float zreion){
    float M0 = 3.0e9;
    float a = 0.17;
    float b = -2.1;
    float d = 2.5;

    float factor = (1.+z)/(1.+zreion);
    float Mcrit = M0 * pow(photHI*1.e12, a) * pow((1.+z)*0.1, b) * pow(1. - factor*factor, d);  
    float fg = pow(2., -Mcrit/Mvir);
    //printf("Mvir = %e\tMcrit = %e\tfg = %f\tzreion = %f\n", Mvir, Mcrit, fg, zreion);
    
    return fg;
}

/*----------------------------------------------------------------------*/
/* RADIATIVE FEEDBACK DESCRIPTION FOLLOWING GNEDIN 2000 & KRAVTSOV 2004 */
/*----------------------------------------------------------------------*/

float fg_radfb_Gnedin(float Mvir, float Temp, float mu, float z, float zion, float zreion, float omega_m)
{
    float fa = get_fa_Gnedin(z, zion, zreion);
    float Mjeans = get_Mjeans(Temp, mu, 0., omega_m);
    
    float sqrt_fa = sqrt(fa);
    float Mcrit = 8. * Mjeans * sqrt_fa*sqrt_fa*sqrt_fa;
    
    float tmp = Mcrit/Mvir;
    float fg = 1./((1. + 0.26*tmp) * (1. + 0.26*tmp) * (1. + 0.26*tmp));
    //printf("z = %f\tMvir = %e\tMcrit = %e\tfg = %f\tzreion = %f\tMjeans = %e\tfa = %f\n", z, Mvir, Mcrit, fg, zreion, Mjeans, fa);

    return fg;
}

/* factor f(a) according to Gnedin 2000 & Kravtsov 2004 */
float get_fa_Gnedin(float z, float zion, float zreion){
    float a = 1./(1.+z);
    float aheat = 0.;
    float aion = 1./(1.+zion);
    float areion = 1./(1.+zreion);
    float alpha = 6.;
    float result = 1.;

    if(a >= areion){
        result = get_fheat(a, aion, aion, aheat, alpha) + get_fion(a, areion, aion) + get_fcool(a, areion);
        //printf("a = %f\t areion = %f\t result = %e %e %e\n", a, areion, get_fheat(a, aion, aion, aheat, alpha), get_fion(a, areion, aion), get_fcool(a, areion));
    }
    else {
        if(a >= aion){
            result = get_fheat(a, aion, aion, aheat, alpha) + get_fion(a, a, aion);
            //printf("a = %f\t aion = %f\t result = %e %e\n", a, areion, get_fheat(a, aion, aion, aheat, alpha), get_fion(a, a, aion));
        }
        else{
            result = get_fheat(a, a, aion, aheat, alpha);
            //printf("a = %f\t a = %f\t result = %e\n", a, areion, result);
            }
        }
    
    return result;
}

float fg_radfb_TempEvol(float Mvir, float Temp, float mu, float z, float zion, float zreion, float omega_m){
    float fa = get_fa(z, zion, zreion, Temp);
    float Mjeans = get_Mjeans(Temp, mu, 0., omega_m);

    float sqrt_fa = sqrt(fa);
    float Mcrit = 8. * Mjeans * sqrt_fa*sqrt_fa*sqrt_fa;

    float tmp = Mcrit/Mvir;
    float fg = 1./((1. + 0.26*tmp) * (1. + 0.26*tmp) * (1. + 0.26*tmp));
    //printf("z = %f\tMvir = %e\tMcrit = %e\tfg = %f\tzreion = %f\tMjeans = %e\tfa = %f\n", z, Mvir, Mcrit, fg, zreion, Mjeans, fa);

    return fg;
}

/* factor f(a) according to own temperature model */
float get_fa(float z, float zion, float zreion, float temp)
{
    float a = 1./(1.+z);
    float arec = 0.000909;
    float adec = 0.00398;
    float aion = 1./(1.+zion);
    float areion = 1./(1.+zreion);
    float result = 1.;

    if(a >= areion){
        result = get_frec(a, adec, arec, temp) + get_fdec(a, aion, adec, temp) + get_fion(a, areion, aion) + get_fcool(a, areion);
        // printf("a = %f\t areion = %f\t result = %e %e %e %e\n", a, areion, get_frec(a, adec, arec), get_fdec(a, aion, adec), get_fion(a, areion, aion), get_fcool(a, areion));
        }
        else{
            if(a >= aion){
                result = get_frec(a, adec, arec, temp) + get_fdec(a, aion, adec, temp) + get_fion(a, a, aion);
                //printf("a = %f\t aion = %f\t result = %e %e %e\n", a, areion, get_frec(a, adec, arec), get_fdec(a, aion, adec), get_fion(a, a, aion));
                }
                else{
                    if(a >= adec){
                        result = get_frec(a, adec, arec, temp) + get_fdec(a, a, adec, temp);
                        //printf("a = %f\t a = %f\t result = %e %e\n", a, areion, get_frec(a, adec, arec), get_fdec(a, a, adec));
                        }
                        else{
                            result = get_frec(a, a, arec, temp);
                        }
                }
            }
    return result;
}

// TODO: it would be nice to add some descriptions of the functions below
float get_frec(float a, float ap, float arec, float temp){
    float TCMB = 2.73;
    float tmp1 = 1. - 2./3. * sqrt(ap/a);
    float tmp2 = 1. - 2./3. * sqrt(arec/a);

    float result = TCMB*3.e-4 / (a * temp) * ( ap * tmp1 - arec * tmp2 );

    return result;
}

float get_fdec(float a, float ap, float adec, float temp){
    float TCMB = 2.73;
    float tmp1 = adec/a;

    float result = TCMB*3.e-4 / temp * tmp1 * (log(ap/adec) - 2.* (sqrt(ap/a) - sqrt(tmp1)));

    return result;
}

float get_fheat(float a, float ap, float aion, float aheat, float alpha){
    float factor1 = 1./(alpha+2.);
    float factor2 = 2./(2*alpha+5.);
    float tmp1 = factor1 - sqrt(ap/a)*factor2;
    float tmp2 = factor1 - sqrt(aheat/a)*factor2;

    float result = 3./a * pow(aion, -alpha) * ( pow(ap, alpha+2)*tmp1 - pow(aheat, alpha+2)*tmp2 );

return result;
}

float get_fion(float a, float ap, float aion){
    float tmp1 = 0.1*ap*ap * (5. - 4.*sqrt(ap/a));
    float tmp2 = 0.1*aion*aion * (5. - 4.*sqrt(aion/a));

    float result = 3./a * (tmp1 - tmp2);

    return result;
}

float get_fcool(float a, float areion){
    float result = areion * (1. - areion/a*(3. - 2.*sqrt(areion/a)));

    return result;
}
