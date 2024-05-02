#ifndef _GLOBAL_H
#define _GLOBAL_H


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define n_data_max 10000
#define n_chi_data_max 10000


#if !defined(GLOBAL)
#  define GLOBAL extern
#endif


typedef char *string;

typedef struct 
{
	string prefix;
	double z; //Evaluation redshift
// Background cosmology:
    double Omm;
    double ns;
    double Omb;
    double h;
    double w;
    double Omc;
    double sigma8;
    double Omnu;   
    double Omw;   
// k table
    string fnamePS;
    string path_Bells;
    double kmin;
    double kmax;
    int Nk;  
      	
	double zbin;  // 

	int chiQuadSteps; // For trapezoidal integration	
	int GLpoints; // For Gaussian Legendre integration	
	int mMax; // Bm moments upto mMax
	int Nell; 
	double ellmin,ellmax;
	int tree_level, writevectors;
	
	int chatty; // =0,1,2.  
	
	int Wg;
	string fWgchi;
	
} cmdline_data, *cmdline_data_ptr;

GLOBAL cmdline_data cmd;


typedef struct
{
	double cpuinit;
	double dx;
    int method_int;
    int quadmethod_int;
	clock_t time;    
	//~ string headline0;
	//~ string headline1;
	//~ string headline2;
	//~ string headline3;
	FILE *outlog;
	char mode[2];
    char fnamePSPath[100];
	double k_data[n_data_max], pkz0_data[n_data_max];
	int n_data;
	double Dpz0;	
	double sigma8; // This is sigma8 for input pk at z=0;	
	//~ int Nk; // number of log spaced k in Bm(k1,k2), Check if used
	//~ double kmin, kmax; // kmin and kmax in Bm(k1,k2). Check if used
	//~ double *kT; // Array of k used for Bm(k1,k2). Check if used	
	//~ double *veckBm0, *veckBm1, *veckBm2;	
	double chi_data[n_chi_data_max], Wg_chi_data[n_chi_data_max];
	int n_chi_data;
	
	 
} global_data, *global_data_ptr;


GLOBAL global_data gd;


typedef struct 
{	
	clock_t time;    
	double Dp, r_sigma, n_eff;
	double z;
	double *chiOfzT, *zT, *DpT;
	double zMin, zMax;
	int Nz;
	double chiBin;  
	double *gLT, *chiforgLT; 
	int NstepsforgL;
} global_vars, *global_vars_ptr;

GLOBAL global_vars gv;


typedef struct 
{		
	double chiMaxInt, chiMinInt;
	double *chiT;
	int chiQuadSteps, Nell;	
	double *chiT_chiint, *zT_chiint, *DpT_chiint, *rsigma_chiint, *neff_chiint;
	double *q_chiint, *kT, *ellT;
	double **BmVectors, **BmVectorsp;
	double ellmin, ellmax;	
} integration_vars, *integration_vars_ptr;

GLOBAL integration_vars iv;




typedef struct 
{	  
	double *r1, *r2, **result;
} zeta_vars, *zeta_vars_ptr;

GLOBAL zeta_vars zv;





#endif  // ! _GLOBAL_H
