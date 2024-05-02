#include "global.h"
#include "procedures.h"

#define eps 0.0001


double Bispec_Takahashi(double k1, double k2, double k3, double z, double Dp, double r_sigma, double n_eff);
double F2(double k1, double k2, double k3, double z, double Dp, double r_sigma);

double F2_tree(double k1, double k2, double k3);
double Bispec_tree(double k1, double k2, double k3, double Dp);
double Bispec_P2(double k1, double k2, double k3, double Dp);
double Bispec_tree_EFT(double k1, double k2, double k3, double Dp, double ctr);


double calcrsigma(double Dp, double kini, double kfin, int Nk);

double sigmaRTH(double r, double kini, double kfin, int Nk);
double sigmaRGaussian(double r, double kini, double kfin, int Nk);
double sigmaRGaussian1stDeriv(double r, double kini, double kfin, int Nk);
double n_eff_func(double r_sigma, double Dp, double kini, double kfin, int Nk);
double sigmam(double r, int j);
double window(double x, int i);

double linear_pkz0(double k); 
double linear_pkz0_data(double k); 
double linear_pkz0_eh(double k);

double Dplusf(double z);
double Dplusf_func(int j, double la, double y[2]);
