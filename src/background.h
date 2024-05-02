#include "global.h"
#include "procedures.h"
#include "functions.h"

void background(double zBin);

void   chiArray_all(double chiOfzT[], double zT[]);
double      HoverH0(double z);
double  chiOfz_func(double z);
double  zOfchi_func(double chi);
double  aOfchi_func(double chi);
double DpOfchi_func(double chi);

double gL(double chi);
double gLDiracDelta(double chi);
double q(double chi);

void chiMaxforInt(void);
void ArraysforChiQuad(void);
void ArraysforChiQuadLog(void);

// for photo-z file
void compute_gL(void);
void read_inputWgchi(void);
double Wg_func(double chi);
double gL_func(double chi);

