#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include <time.h>

#include <fftw3.h>

#include "utils.h"
#include "twobessel.h"


#define m_2PI2 39.4784176043574  // (2 pi)^2

int main(int argc, char const *argv[])
{

	//char inputpath[]="/home/waco/Dropbox/0_3pt/C_code_nbs/making/";
	//char prefix[]="nbs_";
	
	char inputpath[]="../Bell_outputs/";
	char prefix[]="dic24_";
	char pathini[50];
	sprintf(pathini,"%s%s",inputpath,prefix);
	
	char outputpath[]="outputs/";
	char prefixout[]="dic24_";
	char pathout[50]; 
	sprintf(pathout,"%s%s",outputpath,prefixout);	
	


	char file_ells[80];
	char file_Bm0[80];
	char file_Bm1[80];
	char file_Bm2[80];
	char file_Bm3[80];
	
	
	char outfilename_zm0[80], outfilename_zm0_diag[80]; 
	char outfilename_zm1[80], outfilename_zm1_diag[80];
	char outfilename_zm2[80], outfilename_zm2_diag[80]; 
	char outfilename_zm3[80], outfilename_zm3_diag[80];
	
	
	sprintf(file_ells,"%sellArray.txt",pathini);
	sprintf(file_Bm0,"%sBmellsVector_0.txt",pathini);	
	sprintf(file_Bm1,"%sBmellsVector_1.txt",pathini);
	sprintf(file_Bm2,"%sBmellsVector_2.txt",pathini);
	sprintf(file_Bm3,"%sBmellsVector_3.txt",pathini);
	
	

	
	sprintf(outfilename_zm0,"%szeta0.txt",pathout);
	sprintf(outfilename_zm0_diag,"%szeta0_diag.txt",pathout);
	sprintf(outfilename_zm1,"%szeta1.txt",pathout);
	sprintf(outfilename_zm1_diag,"%szeta1_diag.txt",pathout);
	sprintf(outfilename_zm2,"%szeta2.txt",pathout);
	sprintf(outfilename_zm2_diag,"%szeta2_diag.txt",pathout);
	sprintf(outfilename_zm3,"%szeta3.txt",pathout);
	sprintf(outfilename_zm3_diag,"%szeta3_diag.txt",pathout);
	

	
	//~ char filename[] = "../python/Pk_test";
	FILE *INell = fopen(file_ells, "r");

	// double *ell, *fl;
	long Nell = 128;
	double ell[Nell];

	long linenum = 0;
	while(!feof(INell) && (linenum<Nell)) {
		fscanf(INell, "%lg", &ell[linenum]);
		linenum++;
	}
	
	//~ for(long i=0; i<Nell; i++ ){
		//~ printf("%e\n", ell[i]);
	//~ }
	
	double dlnell = log(ell[1]/ell[0]);
	
	

	long size = Nell*Nell;
	long n_data, ndatamax;
    double mat[Nell][Nell];


	config my_config;
	my_config.l1 = 0;
	my_config.l2 = 0;
	my_config.nu1 = 1.01;
	my_config.nu2 = 1.01;
	my_config.c_window_width = 0.25;
	my_config.sys_Flag = 0;


    double Bmell0vector[size];
	FILE *fp;    
	fp = fopen(file_Bm0, "r");    
     
    ndatamax = size;
    
    n_data=0;
;
    
    if (NULL == fp) {
        printf("\n\nfile can't be opened \n\n");
    }
    if(fp!=NULL){   
		while(fscanf(fp, "%lf", &Bmell0vector[n_data])!=EOF){
			n_data++;
			if(n_data>ndatamax) printf("n_data_max should be larger than the number of data lines \n");
		}	  	
	fclose(fp);		
	}
	

 	//~ for(int i=0;i<size; i++){
	//~ printf("%e,  ",Bmell0vector[i]);	
	//~ }	
	

	double **Bmell1ell2;
	Bmell1ell2 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) Bmell1ell2[i] = malloc(Nell * sizeof(double));

	for(int i=0; i<Nell ; i++){
		for(int j=0; j<Nell ; j++){
			Bmell1ell2[i][j] = ell[i]*ell[i]*ell[j]*ell[j]
								*Bmell0vector[i*Nell+j] / m_2PI2;	
		}
	} 


	double *r1, *r2, **result;
	r1 = malloc(Nell * sizeof(double));
	r2 = malloc(Nell * sizeof(double));
	result = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) {
		result[i] = malloc(Nell * sizeof(double));
	}
	
	double smooth_dlnr = dlnell;
	int dimension = 2;
	clock_t start = clock();
	two_Bessel_binave(ell, ell, Bmell1ell2, Nell, Nell, &my_config, smooth_dlnr, dimension, r1, r2, result);
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("time:%f\n", seconds);

	
	FILE *OUT = fopen(outfilename_zm0, "w");	
	for(int i=0; i<Nell; i++) {
		for(int j=0; j<Nell; j++){
	 		fprintf(OUT, "%lg %lg %lg", r1[i], r2[j], result[i][j]);
	 		fprintf(OUT, "\n");
	 	}
	}
	fclose(OUT);

	FILE *OUT2 = fopen(outfilename_zm0_diag, "w");	
	for(int i=0; i<Nell; i++) {
		fprintf(OUT2, "%lg %lg", r1[i], result[i][i]);
		fprintf(OUT2, "\n");
	}
	fclose(OUT2);
	
	free(r1);
	free(r2);
	for (int i = 0; i < Nell; ++i){
		free(Bmell1ell2[i]);
		free(result[i]);
	}
	free(Bmell1ell2);
	free(result);


// Ahora lo mismo para zetam1,zetam2 y zetam3
// Hacer en loop mejor

// zetam1


	config my_config_m1;
	my_config_m1.l1 = 1;
	my_config_m1.l2 = 1;
	my_config_m1.nu1 = 1.01;
	my_config_m1.nu2 = 1.01;
	my_config_m1.c_window_width = 0.25;
	my_config_m1.sys_Flag = 0;

    double Bm1ellvector[size];
	FILE *fpm1;    
	fpm1 = fopen(file_Bm1, "r");    
     
    ndatamax = size;
    
    n_data=0;
    if (NULL == fpm1) {
        printf("\n\nfile can't be opened \n\n");
    }
    if(fpm1!=NULL){   
		while(fscanf(fp, "%lf", &Bm1ellvector[n_data])!=EOF){
			n_data++;
			if(n_data>ndatamax) printf("n_data_max should be larger than the number of data lines \n");
		}	  	
	fclose(fpm1);		
	}

	

	double **Bm1ell1ell2;
	Bm1ell1ell2 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) Bm1ell1ell2[i] = malloc(Nell * sizeof(double));

	for(int i=0; i<Nell ; i++){
		for(int j=0; j<Nell ; j++){
			Bm1ell1ell2[i][j] = ell[i]*ell[i]*ell[j]*ell[j]
								*Bm1ellvector[i*Nell+j] / m_2PI2;	
		}
	} 


	double *r1m1, *r2m1, **resultm1;
	r1m1 = malloc(Nell * sizeof(double));
	r2m1 = malloc(Nell * sizeof(double));
	resultm1 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) {
		resultm1[i] = malloc(Nell * sizeof(double));
	}
	
	//~ double smooth_dlnr = dlnell;
	//~ int dimension = 2;
	start = clock();
	two_Bessel_binave(ell, ell, Bm1ell1ell2, Nell, Nell, &my_config_m1, smooth_dlnr, dimension, r1m1, r2m1, resultm1);
	end = clock();
	seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("time:%f\n", seconds);

	
	FILE *OUTm1 = fopen(outfilename_zm1, "w");	
	for(int i=0; i<Nell; i++) {
		for(int j=0; j<Nell; j++){
	 		fprintf(OUTm1, "%lg %lg %lg", r1m1[i], r2m1[j], resultm1[i][j]);
	 		fprintf(OUTm1, "\n");
	 	}
	}
	fclose(OUTm1);

	FILE *OUT2m1 = fopen(outfilename_zm1_diag, "w");	
	for(int i=0; i<Nell; i++) {
		fprintf(OUT2, "%lg %lg", r1m1[i], resultm1[i][i]);
		fprintf(OUT2, "\n");
	}
	fclose(OUT2m1);
	
	free(r1m1);
	free(r2m1);
	for (int i = 0; i < Nell; ++i){
		free(Bm1ell1ell2[i]);
		free(resultm1[i]);
	}
	free(Bm1ell1ell2);
	free(resultm1);






// zetam2



	config my_config_m2;
	my_config_m2.l1 = 2;
	my_config_m2.l2 = 2;
	my_config_m2.nu1 = 1.01;
	my_config_m2.nu2 = 1.01;
	my_config_m2.c_window_width = 0.25;
	my_config_m2.sys_Flag = 0;


    double Bm2ellvector[size];
	FILE *fpm2;    
	fpm2 = fopen(file_Bm2, "r");    
     
    ndatamax = size;
    
    n_data=0;
    if (NULL == fpm2) {
        printf("\n\nfile can't be opened \n\n");
    }
    if(fpm2!=NULL){   
		while(fscanf(fp, "%lf", &Bm2ellvector[n_data])!=EOF){
			n_data++;
			if(n_data>ndatamax) printf("n_data_max should be larger than the number of data lines \n");
		}	  	
	fclose(fpm2);		
	}

	

	double **Bm2ell1ell2;
	Bm2ell1ell2 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) Bm2ell1ell2[i] = malloc(Nell * sizeof(double));

	for(int i=0; i<Nell ; i++){
		for(int j=0; j<Nell ; j++){
			Bm2ell1ell2[i][j] = ell[i]*ell[i]*ell[j]*ell[j]
								*Bm2ellvector[i*Nell+j] / m_2PI2;	
		}
	} 


	double *r1m2, *r2m2, **resultm2;
	r1m2 = malloc(Nell * sizeof(double));
	r2m2 = malloc(Nell * sizeof(double));
	resultm2 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) {
		resultm2[i] = malloc(Nell * sizeof(double));
	}
	
	//~ double smooth_dlnr = dlnell;
	//~ int dimension = 2;
	start = clock();
	two_Bessel_binave(ell, ell, Bm2ell1ell2, Nell, Nell, &my_config_m2, smooth_dlnr, dimension, r1m2, r2m2, resultm2);
	end = clock();
	seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("time:%f\n", seconds);

	
	FILE *OUTm2 = fopen(outfilename_zm2, "w");	
	for(int i=0; i<Nell; i++) {
		for(int j=0; j<Nell; j++){
	 		fprintf(OUTm2, "%lg %lg %lg", r1m2[i], r2m2[j], resultm2[i][j]);
	 		fprintf(OUTm2, "\n");
	 	}
	}
	fclose(OUTm2);

	FILE *OUT2m2 = fopen(outfilename_zm2_diag, "w");	
	for(int i=0; i<Nell; i++) {
		fprintf(OUT2, "%lg %lg", r1m2[i], resultm2[i][i]);
		fprintf(OUT2, "\n");
	}
	fclose(OUT2m2);
	
	free(r1m2);
	free(r2m2);
	for (int i = 0; i < Nell; ++i){
		free(Bm2ell1ell2[i]);
		free(resultm2[i]);
	}
	free(Bm2ell1ell2);
	free(resultm2);




// zetam3

	config my_config_m3;
	my_config_m3.l1 = 3;
	my_config_m3.l2 = 3;
	my_config_m3.nu1 = 1.01;
	my_config_m3.nu2 = 1.01;
	my_config_m3.c_window_width = 0.25;
	my_config_m3.sys_Flag = 0;

    double Bm3ellvector[size];
	FILE *fpm3;    
	fpm3 = fopen(file_Bm3, "r");    
     
    ndatamax = size;
    
    n_data=0;
    if (NULL == fpm3) {
        printf("\n\nfile can't be opened \n\n");
    }
    if(fpm3!=NULL){   
		while(fscanf(fp, "%lf", &Bm3ellvector[n_data])!=EOF){
			n_data++;
			if(n_data>ndatamax) printf("n_data_max should be larger than the number of data lines \n");
		}	  	
	fclose(fpm3);		
	}

	

	double **Bm3ell1ell2;
	Bm3ell1ell2 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) Bm3ell1ell2[i] = malloc(Nell * sizeof(double));

	for(int i=0; i<Nell ; i++){
		for(int j=0; j<Nell ; j++){
			Bm3ell1ell2[i][j] = ell[i]*ell[i]*ell[j]*ell[j]
								*Bm3ellvector[i*Nell+j] / m_2PI2;	
		}
	} 


	double *r1m3, *r2m3, **resultm3;
	r1m3 = malloc(Nell * sizeof(double));
	r2m3 = malloc(Nell * sizeof(double));
	resultm3 = malloc(Nell * sizeof(double *));
	for(int i=0; i<Nell; i++) {
		resultm3[i] = malloc(Nell * sizeof(double));
	}
	
	//~ double smooth_dlnr = dlnell;
	//~ int dimension = 2;
	start = clock();
	two_Bessel_binave(ell, ell, Bm3ell1ell2, Nell, Nell, &my_config_m3, smooth_dlnr, dimension, r1m3, r2m3, resultm3);
	end = clock();
	seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("time:%f\n", seconds);

	
	FILE *OUTm3 = fopen(outfilename_zm3, "w");	
	for(int i=0; i<Nell; i++) {
		for(int j=0; j<Nell; j++){
	 		fprintf(OUTm3, "%lg %lg %lg", r1m3[i], r2m3[j], resultm3[i][j]);
	 		fprintf(OUTm3, "\n");
	 	}
	}
	fclose(OUTm3);

	FILE *OUT2m3 = fopen(outfilename_zm3_diag, "w");	
	for(int i=0; i<Nell; i++) {
		fprintf(OUT2, "%lg %lg", r1m3[i], resultm3[i][i]);
		fprintf(OUT2, "\n");
	}
	fclose(OUT2m3);
	
	free(r1m3);
	free(r2m3);
	for (int i = 0; i < Nell; ++i){
		free(Bm3ell1ell2[i]);
		free(resultm3[i]);
	}
	free(Bm3ell1ell2);
	free(resultm3);











	return 0;
}
