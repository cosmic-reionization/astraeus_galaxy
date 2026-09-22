#include <cstdint>
//#include <cstddef>

#include "evol_gal.h"
#include "radiative_feedback.h"

#include <Python.h>
#include <numpy/arrayobject.h>

extern "C" {
    static PyObject *astraeus_evolve_gal(PyObject *self, PyObject *args){
    float hubble_h;
    float omega_b;
    float omega_m;
    int32_t localDescID;
    int32_t numProg;
    int32_t snapnumber;
    float Mvir;
    float Mvir_prog;
    float Rvir;
    float Mgas;
    PyArrayObject *igmMetallicity;  // 3-dim array
    float igmDustFraction;
    PyArrayObject *fracMmetalMer;     // 3-dim array
    PyArrayObject *Mmetal;          // 3-dim array
    float Mdust;
    float fracMgasMer;
    float photHI_bg;
    float scalefactor;
    float zreion;
    int32_t reion_flag;
    float temp;
    float mu;
    int32_t radfeedback_model;

    if (!PyArg_ParseTuple(args, "fffiiiffffOfOOfffffiffi", &hubble_h, &omega_b, &omega_m, &localDescID, &numProg, &snapnumber, &Mvir, &Mvir_prog, &Rvir, &Mgas, &igmMetallicity, &igmDustFraction, &fracMmetalMer, &Mmetal, &Mdust, &fracMgasMer, &photHI_bg, &scalefactor, &zreion, &reion_flag, &temp, &mu, &radfeedback_model))
        return NULL;

    // Error checking
    // TODO: consider adding some error checking

    // Get Array data
    float *igmMetallicity_data = (float *)PyArray_DATA(igmMetallicity);
    float *fracMmetalMer_data = (float *)PyArray_DATA(fracMmetalMer);
    float *Mmetal_data = (float *)PyArray_DATA(Mmetal);
    
    evolve_gal(hubble_h, omega_b, omega_m, localDescID, numProg, snapnumber, Mvir, Mvir_prog, Rvir, Mgas, igmMetallicity_data, igmDustFraction, fracMmetalMer_data, Mmetal_data, Mdust, fracMgasMer, photHI_bg, scalefactor, zreion, reion_flag, temp, mu, radfeedback_model);

    Py_RETURN_NONE;
    }

    static PyMethodDef astraeus_methods[] = {
        {"evolve_gal", astraeus_evolve_gal, METH_VARARGS, "Evolve one galaxy."},
        {NULL, NULL, 0, NULL} // Sentinel
    };

    static PyModuleDef astraeus_module = {
        PyModuleDef_HEAD_INIT,                  
        "astraeus",                                     // name of the module (this need to match the name TARGET in the Makefile)
        "Galaxy evolution functions from Astraeus.",    // module documentation
        -1,                                             // size of per-interpreter state of the module, or -1 if the module keeps state in global variables (no return value)
        astraeus_methods
    };

    PyMODINIT_FUNC PyInit_astraeus(void){
        import_array();
        return PyModule_Create(&astraeus_module);   // import module
    }
}
