#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "evol_gal.h"
#include "radiative_feedback.h"
//#include "metallicity.h"
#include "starformation_SNfeedback.h"


// -------- global variables --------
// TODO: check with Anne why this is needed and if it's indeed (1*u.yr).to('s') 
const float YrInSec = 31557600.;

/*-----------------------------------------------------------------------*/
/* This function evolves the properties of a galaxy for one time step    */
/* Version modified by Michele Bianco (11.09.2026) for wrapping for a Python interface */
/*-----------------------------------------------------------------------*/
//void do_evolve_gal()


void evolve_gal(const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &localDescID, const int32_t &numProg, const int32_t &snapnumber, float &Mvir, float &Mvir_prog, float &Rvir, float &Mgas, float *igmMetallicity, float &igmDustFraction, float *fracMmetalMer, float *Mmetal, float &Mdust, float &fracMgasMer, const float &photHI_bg, const float &scalefactor, const float &zreion, const int32_t &reion_flag, const float &temp, const float &mu, const int32_t &radfeedback_model){

    /* adjust DM mass to continuously rise */
    correct_Mvir_to_continously_rise(Mvir, Mvir_prog, Rvir);
    //printf("Mvir = %e\t Rvir = %e\t Mvir_prog = %e\n\n", Mvir, Rvir, Mvir_prog);

    /* get metals from inflows & progenitors */
    //get_accreted_and_merged_metals(omega_b, omega_m, numProg, snapnumber, igmMetallicity, igmDustFraction, Mvir, Mvir_prog, fracMmetalMer, Mmetal, Mdust);
    //printf("Mmetal = %e %e %e\nMdust = %e\nfracMmetalMer = %e %e %e\n\n", Mmetal[0], Mmetal[1], Mmetal[2], Mdust, fracMmetalMer[0], fracMmetalMer[1], fracMmetalMer[2]);

    /* get gas from inflows & progenitors */
    //get_accreted_and_merged_gas(omega_b, omega_m, numProg, Mgas, Mvir, Mvir_prog, fracMgasMer);
    //printf("Mgas = %e\nfracMgasMer = %e\n\n", Mgas, fracMgasMer);

    /* get metals and gas from inflows & progenitors */
    get_accreted_and_merged(omega_b, omega_m, numProg, snapnumber, igmMetallicity, igmDustFraction, Mvir, Mvir_prog, Mgas, fracMgasMer, Mmetal, fracMmetalMer, Mdust);

    //TODO: Here I merged the apply_radFeedback_and_update_zreion and apply_radFeedback_to_metals functions into one, since both shared common inputs
    /* radiative feedback (how much gas can the galaxy keep) */
    //apply_radFeedback_and_update_zreion(hubble_h, omega_b, omega_m, snapnumber, Mvir, Mgas, photHI_bg, scalefactor, zreion, reion_flag, temp, mu, radfeedback_model);

    //if(simParam->metals == 1)
    //apply_radFeedback_to_metals(omega_b, omega_m, Mvir, Mgas, fg, Mmetal, Mdust);

    apply_radFeedback_and_metals(hubble_h, omega_b, omega_m, snapnumber, Mvir, Mgas, Mmetal, Mdust, photHI_bg, scalefactor, zreion, reion_flag, temp, mu, radfeedback_model);
    
    /* saving initial gas mass (before SN feedback) */
    float MgasIni = Mgas;

    /* saving initial metal masses (before SN feedback) */
    float MmetalIni[3];
    MmetalIni[0] = Mmetal[0];
    MmetalIni[1] = Mmetal[1];
    MmetalIni[2] = Mmetal[2];

    /* star formation & SN feedback */
    //do_starFormation_and_SNfeedback(tf, ti, dt_deltaTimeInMyr, Mvir, Rvir, Mgas, fracMgasMer, MgasIni, Mstar, hubble_h, omega_b, omega_m, snapnumber, scalefactor, temp, mu, zreion, reion_flag, radfeedback_model, delayedSNfeedback, FS, FW, SNenergy);

    //produce_and_eject_metals(simParam, thisGal, descGal);

    //TODO: return everything that need to be used in Python

}

/*---------------------------------------------------------------------*/
/* Correct Mvir & Rvir if Mvir is not continuously rising in trees     */
/*---------------------------------------------------------------------*/
void correct_Mvir_to_continously_rise(float &Mvir, float &Mvir_prog, float &Rvir){
    if(Mvir < Mvir_prog){
        Rvir *= pow(Mvir_prog/Mvir, 0.3333);
        Mvir = Mvir_prog;
    } 
}

/*---------------------------------------------------------------------*/
/* Get accreted & merged metal mass                                    */
/*---------------------------------------------------------------------*/
void get_accreted_and_merged(const float &omega_b, const float &omega_m, const int32_t &numProg, const int32_t &snapnumber, float *igmMetallicity, float &igmDustFraction, float &Mvir, float &Mvir_prog, float &Mgas, float &fracMgasMer, float *Mmetal, float *fracMmetalMer, float &Mdust){
    float MgasAcc = 0.;
    float MmetalAcc = 0., MmetalAccO = 0., MmetalAccFe = 0.;
    float dustfractionIgm = 0.;

    // Compute average IGM metallicity from previous snap (metallicity defined in Python)
    float metallicityIgm   = igmMetallicity[0];
    float metallicityIgmO  = igmMetallicity[1];
    float metallicityIgmFe = igmMetallicity[2];
    dustfractionIgm = igmDustFraction;
    printf("before %e\n", Mgas);

    // Get gas from INFLOWS and PROGENITORS
    if(numProg == 0){
        Mgas = Mvir * omega_b / omega_m;
        fracMgasMer = 0.;

        MmetalAcc = metallicityIgm * Mgas;
        MmetalAccO = metallicityIgmO * Mgas;
        MmetalAccFe = metallicityIgmFe * Mgas;
        Mmetal[0] = MmetalAcc;
        Mmetal[1] = MmetalAccO;
        Mmetal[2] = MmetalAccFe;
        Mdust  = dustfractionIgm * Mgas;

    }
    else if(Mvir_prog < Mvir){
        MgasAcc = (Mvir - Mvir_prog) * omega_b / omega_m;
        Mgas += MgasAcc;
        fracMgasMer = 1.0 - MgasAcc/Mgas;

        MmetalAcc = metallicityIgm * MgasAcc;
        MmetalAccO = metallicityIgmO * MgasAcc;
        MmetalAccFe = metallicityIgmFe * MgasAcc;
        Mmetal[0] += MmetalAcc;
        Mmetal[1] += MmetalAccO;
        Mmetal[2] += MmetalAccFe;
        Mdust  += dustfractionIgm * MgasAcc;
    }

    fracMmetalMer[0] = 1. - (MmetalAcc / Mmetal[0]);
    fracMmetalMer[1] = 1. - (MmetalAccO / Mmetal[1]);
    fracMmetalMer[2] = 1. - (MmetalAccFe / Mmetal[2]);
    printf("after %e\n", Mgas);
}


void get_accreted_and_merged_metals(const float &omega_b, const float &omega_m, const int32_t &numProg, const int32_t &snapnumber, float *igmMetallicity, float &igmDustFraction, float &Mvir, float &Mvir_prog, float *fracMmetalMer, float *Mmetal, float &Mdust){
    float MgasAcc = 0.;
    float MmetalAcc = 0., MmetalAccO = 0., MmetalAccFe = 0.;
    float metallicityIgm = 0., metallicityIgmO = 0., metallicityIgmFe = 0.;
    float dustfractionIgm = 0.;

    // Compute average IGM metallicity from previous snap 
    if(snapnumber != 0){
        metallicityIgm   = igmMetallicity[0]; //simParam->meanMetallicityIgm[thisGal->snapnumber-1];
        metallicityIgmO  = igmMetallicity[1]; //simParam->meanMetallicityIgmO[thisGal->snapnumber-1];
        metallicityIgmFe = igmMetallicity[2]; //simParam->meanMetallicityIgmFe[thisGal->snapnumber-1];
        dustfractionIgm = igmDustFraction; //simParam->meanDustfractionIgm[thisGal->snapnumber-1];
    }

    // Get gas from INFLOWS and PROGENITORS
    if(numProg == 0){
        MgasAcc = Mvir * omega_b / omega_m;
        MmetalAcc = metallicityIgm * MgasAcc;
        MmetalAccO = metallicityIgmO * MgasAcc;
        MmetalAccFe = metallicityIgmFe * MgasAcc;
        Mmetal[0] = MmetalAcc;
        Mmetal[1] = MmetalAccO;
        Mmetal[2] = MmetalAccFe;
        Mdust  = dustfractionIgm * MgasAcc;
    }
    else if(Mvir_prog < Mvir){
        MgasAcc = (Mvir - Mvir_prog) * omega_b / omega_m;
        MmetalAcc = metallicityIgm * MgasAcc;
        MmetalAccO = metallicityIgmO * MgasAcc;
        MmetalAccFe = metallicityIgmFe * MgasAcc;
        Mmetal[0] += MmetalAcc;
        Mmetal[1] += MmetalAccO;
        Mmetal[2] += MmetalAccFe;
        Mdust  += dustfractionIgm * MgasAcc;
    }

    fracMmetalMer[0] = 1. - (MmetalAcc / Mmetal[0]);
    fracMmetalMer[1] = 1. - (MmetalAccO / Mmetal[1]);
    fracMmetalMer[2] = 1. - (MmetalAccFe / Mmetal[2]);
}

/*---------------------------------------------------------------------*/
/* Get accreted & merged gas mass & fraction of merged gas             */
/*---------------------------------------------------------------------*/
void get_accreted_and_merged_gas(const float &omega_b, const float &omega_m, const int32_t &numProg, float &Mgas, float &Mvir, float &Mvir_prog, float &fracMgasMer){
    // Gas mass before the merge
    //float MgasMer = 0.;
    float MgasAcc = 0;

    /* Equation (1) and (2) in https://arxiv.org/pdf/2004.08401 */
    if(numProg == 0){
        // TODO: there is several repetition of this formula in the other functions. Maybe we can instead pass the Mgas directly
        Mgas = Mvir * omega_b / omega_m;
        fracMgasMer = 0.;
    }
    else if(Mvir_prog < Mvir){
        MgasAcc = (Mvir - Mvir_prog) * omega_b / omega_m;
        Mgas += MgasAcc;
        fracMgasMer = 1.0 - MgasAcc/Mgas;
        //MgasMer = Mgas;
        //Mgas += (Mvir - Mvir_prog) * omega_b / omega_m;
        //fracMgasMer = MgasMer / Mgas;
    }
}

/*---------------------------------------------------------------------*/
/* Compute reduction of gas mass for star formation due to radiative   */
/* feedback and track when a galaxy's environment was reionized        */
/*---------------------------------------------------------------------*/
/*
void apply_radFeedback_and_update_zreion(const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &snapnumber, float &Mvir, float &Mgas, const float &photHI_bg, const float &scalefactor, const float &zreion, const int32_t &reion_flag, const float &temp, const float &mu, const int32_t &radfeedback_model){
    float fg = 1.;
    float MgasMax = 0.;

    //if(radfeedback == 1){
    // TODO: Should we always assume radiative feedback is on? since we are using in pyc2ray
    if(true){
        fg = fmin(1., calc_radiative_feedback(hubble_h, omega_m, snapnumber, Mvir, photHI_bg, scalefactor, zreion, reion_flag, temp, mu, radfeedback_model));
        
        assert(fg <= 1.);
        
        MgasMax = fg * omega_b / omega_m * Mvir;
        if(Mgas > MgasMax){
            Mgas = MgasMax;
        } //TODO: question for Anna. Here we are reducing the gas mass of the galaxy if it is above the maximum allowed by radiative feedback. But what if below?
    }
    /* TODO: ask Anne what does this means and if it's needed. From what I remember from our conversation this is needed by the reionisation of astraeus but if we embended in pyc2ray is not needed
    else if(reion_flag == 1){
        pass_zreion_to_descGal(thisGal, descGal, simParam);
    }
}


void apply_radFeedback_to_metals(const float &omega_b, const float &omega_m, float &Mvir, float &Mgas, float &fg, float *Mmetal, float &Mdust)
{
    float metallicityIni = 0., metallicityIniO = 0., metallicityIniFe = 0.;
    float dustfractionIni = 0.;
    float MgasMax = 0.;

    //if(simParam->radfeedback == 1){}
    // TODO: Should we always assume radiative feedback is on? since we are using in pyc2ray

    if(Mgas > 0.){
        metallicityIni    = Mmetal[0] / Mgas;
        metallicityIniO  = Mmetal[1] / Mgas;
        metallicityIniFe = Mmetal[2] / Mgas;
        dustfractionIni    = Mdust  / Mgas; 
    }

    // TODO: there is a some repetition of this formula in apply_radFeedback_and_update_zreion. Maybe we should consider merge the two.
    MgasMax = fg * omega_b / omega_m * Mvir;
    
    if(Mgas > MgasMax){
        Mmetal[0] = metallicityIni * MgasMax;
        Mmetal[1] = metallicityIniO * MgasMax;
        Mmetal[2] = metallicityIniFe * MgasMax; 
        Mdust  = dustfractionIni * MgasMax;
    }
    
}
*/

void apply_radFeedback_and_metals(const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &snapnumber, float &Mvir, float &Mgas, float *Mmetal, float &Mdust, const float &photHI_bg, const float &scalefactor, const float &zreion, const int32_t &reion_flag, const float &temp, const float &mu, const int32_t &radfeedback_model){
    // TODO: this function is a merge of the apply_radFeedback_and_update_zreion and apply_radFeedback_to_metals functions. The reason is that both functions share common inputs and outputs, so it makes sense to combine them into one function to avoid code duplication and improve efficiency.
    // TODO: look at the original functions (here above) and discuss with Anne.
    
    // variables from apply_radFeedback_and_update_zreion
    float fg = 1.;
    float MgasMax = 0.;

    // variables from apply_radFeedback_to_metals
    float metallicityIni = 0., metallicityIniO = 0., metallicityIniFe = 0.;
    float dustfractionIni = 0.;

    fg = fmin(1., calc_radiative_feedback(hubble_h, omega_m, snapnumber, Mvir, photHI_bg, scalefactor, zreion, reion_flag, temp, mu, radfeedback_model));
    printf("fg = %e\n", fg);

    assert(fg <= 1.);
    
    MgasMax = fg * omega_b / omega_m * Mvir;
    
    // Eqution (3) in paper I
    Mgas = fminf(Mgas, MgasMax);

    // If you have feedback you might reduce the amought of gas, but you want to keep the same fraction metallicity
    if(Mgas > 0.){
        metallicityIni = Mmetal[0] / Mgas;
        metallicityIniO = Mmetal[1] / Mgas;
        metallicityIniFe = Mmetal[2] / Mgas;
        dustfractionIni = Mdust  / Mgas; 
    }
    
    //TODO: because of the fminf before Mgas is smaller or equal MgasMax, but never larger. So this condition is never met and Mmetal vector is not updated. Eventually we can delete these two if conditions.
    if(Mgas > MgasMax){
        Mmetal[0] = metallicityIni * MgasMax;
        Mmetal[1] = metallicityIniO * MgasMax;
        Mmetal[2] = metallicityIniFe * MgasMax; 
        Mdust  = dustfractionIni * MgasMax;        
    }
}

/*---------------------------------------------------------------------*/
/* Compute stellar mass formed, track stellar mass history and gas     */
/* ejected due to SN feedback                                          */
/*---------------------------------------------------------------------*/
void do_starFormation_and_SNfeedback(const float &tf, const float &ti, const float &dt_deltaTimeInMyr, float &Mvir, float &Rvir,float &Mgas, float &fracMgasMer, float &MgasIni, float &Mstar, const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &snapnumber, const float &scalefactor, const float &temp, const float &mu, const float &zreion, const int32_t &reion_flag, const int32_t &radfeedback_model, const int &delayedSNfeedback, float &FS, float &FW, float *SNenergy){
    float fej = 1.;
    float corrFactorTimeStep = get_corrFactorTimeStep(tf, ti, dt_deltaTimeInMyr, snapnumber);
    float fstar = FS * corrFactorTimeStep;
    float feff = fstar;
    float newMstar = 0.;

    /* supernova feedback */
    fej = calc_SNejection_efficiency(delayedSNfeedback, FW, SNenergy, snapnumber, Mvir, Rvir, scalefactor, fracMgasMer, MgasIni, corrFactorTimeStep);
    
    /* star formation */
    feff = fmin(fstar, fej);
    newMstar = Mgas*feff;
    
    /* total stellar mass and final gas mass */
    Mstar += newMstar;
    // To be sure nothing unphysical appens
    if(fej > 0.){
        Mgas = (Mgas - newMstar) * calc_remaining_gasfraction(delayedSNfeedback, FW, SNenergy, snapnumber, Mvir, Rvir, scalefactor, fracMgasMer, MgasIni, corrFactorTimeStep, fej, feff);
    }
    else{
        Mgas = 0.;
    }

    /* storing stellar mass history of the galaxy */
    //TODO: return Mstar, Mgas, newMstar
}

float get_corrFactorTimeStep(const float tf, const float ti, const float dt_deltaTimeInMyr, int snap){
    float corrFactorTimeStep = 1.;
    if(snap > 0){
        // TODO: delta t of the snapshots
        // TODO: numerical value should be years in sec
        corrFactorTimeStep = (tf - ti) * 3.170979e-14 / dt_deltaTimeInMyr;
    }
    return corrFactorTimeStep;
}

void produce_and_eject_metals(const int snapnumber, float MgasIni, float *Mmetal, float &Mvir, float &Mgas, float &Mdust, float &newMstar, float *MmetalEj, float &MdustEj){
    //int snap = thisGal->snapnumber;
    float metallicityIni = 0., metallicityIniO = 0., metallicityIniFe = 0., dustfractionIni = 0.;
    float metallicity = 0., metallicityO = 0., metallicityFe = 0., dustfraction = 0.;
    //TODO: ask Anne to clarify the function stellarmasshistory
    //float newMstar = thisGal->stellarmasshistory[snap];
    float sfr = 0.;
    float tmpMgas = MgasIni, MgasEj = 0.;

    float G_t       = 0.;    /* (G(t) x Deltat) in Ucci et al. (2021) */
    float e_Z       = 0.;    /* (e_Z  x Deltat) in Ucci et al. (2021) */
    float e_Z_O     = 0.;    /* (e_Z  x Deltat) in Ucci et al. (2021) for OXYGEN */
    float e_Z_Fe    = 0.;    /* (e_Z  x Deltat) in Ucci et al. (2021) for IRON */
    float y_d       = 0.;    /* (y_d  x Deltat) in Ucci et al. (2021) for dust */
    float zetaej    = 0.;    /* metal loading factor for ejection */
    
    /* Compute initial metallicities */
    if(MgasIni != 0){
        metallicityIni   = Mmetal[0] / MgasIni;
        metallicityIniO  = Mmetal[1] / MgasIni;
        metallicityIniFe = Mmetal[2] / MgasIni;
        dustfractionIni  = Mdust / MgasIni;
    }
    
    //TODO: ask Anne to clarify the function stellarmasshistory
    //start_metalmasshistory(thisGal);
    
    /* Compute SFR in galaxy [ M_sun yr^-1 ]  */
    if(snapnumber > 0)
        // TODO: invdeltat is the snapshots 1/dt=(tf-ti)^-1
        sfr = YrInSec * newMstar; //* simParam->invdeltat[snapnumber];

    /* Compute metal and gas masses from SN explosions */
    tmpMgas -= newMstar; 
    Mmetal[0] -= metallicityIni   * newMstar;
    Mmetal[1] -= metallicityIniO  * newMstar;
    Mmetal[2] -= metallicityIniFe * newMstar;
    Mdust     -= dustfractionIni  * newMstar;

    //TODO: isolate this function (see metallicity.cpp)
    //compute_metal_mass(thisGal, simParam, sfr, (1. / simParam->invdeltat[snapnumber]), &G_t, &e_Z, &e_Z_O, &e_Z_Fe, &y_d);
    
    /* fix when gsl integration or so fails, returning NAN */
    if(snapnumber < 1){
        if(std::isnan(G_t)) G_t = 0.;
        if(std::isnan(e_Z)) e_Z = 0.;
        if(std::isnan(e_Z_O)) e_Z_O = 0.;
        if(std::isnan(e_Z_Fe)) e_Z_Fe = 0.;
        if(std::isnan(y_d)) y_d = 0.;
    }

    tmpMgas += G_t + e_Z;
    Mmetal[0] += e_Z;
    Mmetal[1] += e_Z_O;
    Mmetal[2] += e_Z_Fe;
    Mdust     += y_d;

    // Assign metal and gass masses released to outputs
    //TODO: this need to be returned because you use it for (eventually) totall luminosity history for the SED metallicity evolution
    /*
    MgasNew = G_t + e_Z;
    MmetalNew[0] = e_Z;
    MmetalNew[1] = e_Z_O;
    MmetalNew[2] = e_Z_Fe; 
    */
    
    /* get metallicities before ejection */
    metallicity   = Mmetal[0] / tmpMgas;
    metallicityO  = Mmetal[1] / tmpMgas;
    metallicityFe = Mmetal[2] / tmpMgas;
    dustfraction  = Mdust / tmpMgas;
    
    
    /* the ejected masses */
    float ejectedGasFraction = 1. - Mgas / (MgasIni - newMstar);
    if(MgasIni - newMstar == 0.){
        ejectedGasFraction = 0.;
    }
    assert(ejectedGasFraction <= 1. && ejectedGasFraction >= 0.);
    
    if(Mgas > 0.){
        // Note that tmpMgas now is (thisGal->MgasIni - newMstar + thisGal->MgasNew)
        MgasEj = tmpMgas * ejectedGasFraction;
        // TODO: return this vector + dust
        MmetalEj[0]   = MgasEj * metallicity;
        MmetalEj[1]  = MgasEj * metallicityO;
        MmetalEj[2] = MgasEj * metallicityFe;
        MdustEj    = MgasEj * dustfraction;

        // Metal loading factor for ejection
        //TODO: for future idea on more efficient metals ejection
        zetaej = 1.0;
        //zetaej = simParam->metal_ejectLoadingFactor;
        //zetaej = 2 + pow((8 / log10(Mvir)),5);
        //zetaej = 1.0;//0.5 + pow((1e7 / Mvir),0.3);
        
        tmpMgas -= MgasEj;
        Mmetal[0] -= (zetaej * MmetalEj[0]);
        Mmetal[1] -= (zetaej * MmetalEj[1]);
        Mmetal[2] -= (zetaej * MmetalEj[2]);
        Mdust     -= (zetaej * MdustEj);
    }
    else{
        MgasEj = tmpMgas;
        MmetalEj[0]   = Mmetal[0];
        MmetalEj[1]  = Mmetal[1];
        MmetalEj[2] = Mmetal[2];
        MdustEj    = Mdust;

        tmpMgas = 0.;
        Mmetal[0] = 0.;
        Mmetal[1] = 0.;
        Mmetal[2] = 0.; 
        Mdust  = 0.;
    }

    // can this be changed to a maximum function
    Mmetal[0] = fmaxf(Mmetal[0], 0);
    Mmetal[1] = fmaxf(Mmetal[1], 0);
    Mmetal[2] = fmaxf(Mmetal[2], 0);
    Mdust = fmaxf(Mdust, 0);

    // Save metal ejected to output
    // TODO: return MgasEj, MmetalEj and MdustEj
    
    //TODO: remove tmpMgas and instead add/substract diretcly
    Mgas = tmpMgas;

    /* storing metal mass history of the galaxy */
    // TODO: return MgasNew, MmetalNew (see above)
}