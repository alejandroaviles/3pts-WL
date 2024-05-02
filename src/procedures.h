/*==============================================================================
 HEADER: procedures.h				[gsm]
 ==============================================================================*/

#ifndef _procedures_h
#define _procedures_h

// Main procedures //
void InputParams(void);
void Initial(void);
void MainLoop(void);
void Tests(void);
//~ void StartRun(string, string, string, string);
//~ void StartOutput(void);
void EndRun(void);



// integration
void allocate_iv(void);
void Bmell(void);
void Bm(double chi, double z, double Dp, double r_sigma, double n_eff);

// others
void BmKspace(int Maxm, double kmin, double kmax, int Nk, int GLpoints, double z, double Dp, double r_sigma, double n_eff);


// routines
void read_inputpk(void); 
void gaussleg(double x1, double x2, double x[], double w[], int n);
//~ void GLvectors(void);
double interpolationlog(double x, double xT[], double yT[], int n_data);
double interpolation1(double x, double xT[], double yT[], int n_data);


void write(void);

void free_variables(void);


// Tests:  in tests.c
void tests2(void);
void vectortomatrix(void);
void testTakahashiBispectrum(void);






#endif // ! _procedures_h
