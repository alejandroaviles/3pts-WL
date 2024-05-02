#include "zetam.h"

#include <math.h>
#include <complex.h>
#include <string.h>

#include <time.h>
#include <fftw3.h>

void allocate_zv(void)
{
	//~ double *r1, *r2, **result;
	zv.r1 = malloc(iv.Nell * sizeof(double));
	zv.r2 = malloc(iv.Nell * sizeof(double));
	zv.result = malloc(iv.Nell * sizeof(double *));
	for(int i=0; i<iv.Nell; i++) {
		zv.result[i] = malloc(iv.Nell * sizeof(double));
	} 
	
}



void get_zetam()
{


	char outfilename_zm0[80]; 
	char outfilename_zm1[80];
	char outfilename_zm2[80]; 
	char outfilename_zm3[80];
	
	char outfilename_thetavector[80];
	char outfilename_zm0_matrix[80];	

	char pathout[50]; 
	sprintf(pathout,"%s/%s",cmd.path_Bells,cmd.prefix);
	
	

	sprintf(outfilename_zm0,"%szeta0.txt",pathout);
	sprintf(outfilename_zm1,"%szeta1.txt",pathout);
	sprintf(outfilename_zm2,"%szeta2.txt",pathout);
	sprintf(outfilename_zm3,"%szeta3.txt",pathout);

	sprintf(outfilename_zm0_matrix,"%s/%szetam0_matrix.txt",cmd.path_Bells,cmd.prefix);
	sprintf(outfilename_thetavector,"%s/%stheta_array.txt",cmd.path_Bells,cmd.prefix);

	printf("\n\n %s \n\n",outfilename_zm0);



	//~ for(int i=0; i<iv.Nell ; i++){		
		//~ printf("%15e\n", iv.ellT[i]);		
	//~ }

	//~ for(int i=0; i<200 ; i++){		
		//~ printf("%e\n", iv.BmVectorsp[0][i]);		
	//~ }


	int multipole;

	double dlnell = log(iv.ellT[1]/iv.ellT[0]);
	
	int Nell=iv.Nell;
	

	long size = Nell*Nell;
	long n_data, ndatamax;
    double mat[Nell][Nell];

	
	multipole=0;
	
	config config_m0;
		config_m0.l1 = 0;
		config_m0.l2 = 0;
		config_m0.nu1 = 1.01;
		config_m0.nu2 = 1.01;
		config_m0.c_window_width = 0.25;
		config_m0.sys_Flag = 0;
	

	
	double **Bm0ell1ell2;
	Bm0ell1ell2 = malloc(iv.Nell * sizeof(double *));
	for(int i=0; i<iv.Nell; i++) Bm0ell1ell2[i] = malloc(iv.Nell * sizeof(double));

	for(int i=0; i<iv.Nell ; i++){
		for(int j=0; j<iv.Nell ; j++){
			Bm0ell1ell2[i][j] = iv.ellT[i]*iv.ellT[i]*iv.ellT[j]*iv.ellT[j]
								* iv.BmVectorsp[multipole][i*iv.Nell+j] / m_2PI2;	
		}
	} 
	


	
	double smooth_dlnr = dlnell;
	int dimension = 2;
	clock_t start = clock();
	two_Bessel_binave(iv.ellT, iv.ellT, Bm0ell1ell2, iv.Nell, iv.Nell, &config_m0, smooth_dlnr, dimension, zv.r1, zv.r2, zv.result);
	//~ printf("%lg %lg %lg", zv.r1[0], zv.r2[0], zv.result[0][0]);
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("time zetam_0: %f \n", seconds);

	FILE *OUTtheta = fopen(outfilename_thetavector, "w");	
	for(int i=0; i<iv.Nell; i++) fprintf(OUTtheta, "%lg\n", zv.r1[i]);
	fclose(OUTtheta);
	

	
	FILE *OUT = fopen(outfilename_zm0, "w");	
	for(int i=0; i<iv.Nell; i++) {
		for(int j=0; j<iv.Nell; j++){
	 		fprintf(OUT, "%lg %lg %lg", zv.r1[i], zv.r2[j], zv.result[i][j]);
	 		fprintf(OUT, "\n");
	 	}
	}
	fclose(OUT);


	FILE *OUTm0 = fopen(outfilename_zm0_matrix, "w");	
	for(int i=0; i<iv.Nell; i++) {
		for(int j=0; j<iv.Nell; j++){
	 		fprintf(OUTm0, "%lg  ", zv.result[i][j]);
	 	}
	 	if(i!=iv.Nell-1) fprintf(OUTm0, "\n");
	}
	fclose(OUTm0);



	free(zv.r1);
	free(zv.r2);
	for (int i = 0; i < Nell; ++i){
		free(Bm0ell1ell2[i]);
		free(zv.result[i]);
	}
	free(Bm0ell1ell2);
	free(zv.result);



}












/////////////////////////////////////////////
/* free variables */
/////////////////////////////////////////////

	void free_variables(void)
	{	
			for (int m=0; m<cmd.mMax+1; m++){
				free(iv.BmVectorsp[m]);
				free(iv.BmVectors [m]);
			}	
			free(iv.BmVectorsp);
			free(iv.BmVectors );	
	}











