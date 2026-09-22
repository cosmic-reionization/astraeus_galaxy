import numpy as np
import astraeus

from astropy.cosmology import Planck18 as cosmo

igm_metallicity = np.array([1e-3, 1e-3, 1e-3], dtype=np.float32)
frac_mmetal_mer = np.zeros(3, dtype=np.float32)
mmetal = np.array([1e5, 1e5, 1e5], dtype=np.float32)

Mvir, Mvir_prog = float(1e10), float(9e9)
Mgas = Mvir * cosmo.Om0 / cosmo.Ob0 * 1e-3

print("Before:")
print(" Mvir = %e" %Mvir)
print(" Mvir_prog = %e" %Mvir_prog)
print(" Mgas = %e" %Mgas)
print(" Metals = ", mmetal)

print("\n...some internal printf:")
#astraues.do_evolve_gal(halo[i: 1+100])
astraeus.evolve_gal(cosmo.h, cosmo.Ob0, cosmo.Om0, 
                    1,               # localDescID
                    2,               # numProg
                    10,              # snapnumber
                    Mvir,          # Mvir
                    Mvir_prog,           # Mvir_prog
                    100.0,           # Rvir
                    Mgas,           # Mgas
                    igm_metallicity, # igmMetallicity
                    0.1,             # igmDustFraction
                    frac_mmetal_mer, # fracMmetalMer
                    mmetal,          # Mmetal
                    1.0e7,           # Mdust
                    0.1,              # fracMgasMer
                    1e-12,            # photHI_bg
                    1.0,              # scalefactor
                    8.0,              # zreion    
                    1,                # reion_flag
                    1e4,              # temp
                    0.6,              # mu
                    1                 # radfeedback_model (1: Sobacchi & Mesinger 2013, 2: Gnedin 2000 & Kravtsov 2004, 3: fg_radfb_Mjeans 4: fg_radfb_Mcool, 5: fg_radfb_TempEvol)
                    )

print("\nAfter:")
print(" Mvir = %e" %Mvir)
print(" Mgas = %e" %Mgas)
print(" Metals = ", mmetal)
