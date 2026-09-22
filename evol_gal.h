#pragma once
#include <cstdint>

void evolve_gal(const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &localDescID, const int32_t &numProg, const int32_t &snapnumber, float &Mvir, float &Mvir_prog, float &Rvir, float &Mgas, float *igmMetallicity, float &igmDustFraction, float *fracMmetalMer, float *Mmetal, float &Mdust, float &fracMgasMer, const float &photHI_bg, const float &scalefactor, const float &zreion, const int32_t &reion_flag, const float &temp, const float &mu, const int32_t &radfeedback_model);

void correct_Mvir_to_continously_rise(float &Mvir, float &Mvir_prog, float &Rvir);

void get_accreted_and_merged_metals(const float &omega_b, const float &omega_m, const int32_t &numProg, const int32_t &snapnumber, float *igmMetallicity, float &igmDustFraction, float &Mvir, float &Mvir_prog, float *fracMmetalMer, float *Mmetal, float &Mdust);

void get_accreted_and_merged_gas(const float &omega_b, const float &omega_m, const int32_t &numProg, float &Mgas, float &Mvir, float &Mvir_prog, float &fracMgasMer);

void get_accreted_and_merged(const float &omega_b, const float &omega_m, const int32_t &numProg, const int32_t &snapnumber, float *igmMetallicity, float &igmDustFraction, float &Mvir, float &Mvir_prog, float &Mgas, float &fracMgasMer, float *Mmetal, float *fracMmetalMer, float &Mdust);

//void apply_radFeedback_and_update_zreion(const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &snapnumber, float &Mvir, float &Mgas, const float &photHI_bg, const float &scalefactor, const float &zreion, const int32_t &reion_flag, const float &temp, const float &mu, const int32_t &radfeedback_model);

//void apply_radFeedback_to_metals(const float &omega_b, const float &omega_m, float &Mvir, float &Mgas, float &fg, float *Mmetal, float &Mdust);

void apply_radFeedback_and_metals(const float &hubble_h, const float &omega_b, const float &omega_m, const int32_t &snapnumber, float &Mvir, float &Mgas, float *Mmetal, float &Mdust, const float &photHI_bg, const float &scalefactor, const float &zreion, const int32_t &reion_flag, const float &temp, const float &mu, const int32_t &radfeedback_model);

float get_corrFactorTimeStep(const float tf, const float ti, const float dt_deltaTimeInMyr, int snap);
