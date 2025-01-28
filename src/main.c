#define GLOBAL

#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "procedures.h"
#include "functions.h"
#include "background.h"
#include "zetam.h"

int main(int argc, string argv[])
{
	clock_t start_t, end_t;
	start_t = clock();
	//~ InitParam(argv, defv);

	InputParams();
	Initial();
	MainLoop();
	//~ Tests();
	free_variables();
	//~ EndRun();
	end_t = clock(); 
	printf("Total time: %f\n", (double)(end_t - start_t) / CLOCKS_PER_SEC  );
	return 0;
}



void InputParams()
{

	cmd.z= 1.0334;   //No hace nada    
	// cosmological parameters WMAP 2009: 
	cmd.h=0.7;     // Hubble parameter 
	cmd.sigma8=0.82; // sigma 8
	cmd.Omb=0.046;   // Omega baryon
	cmd.Omc=0.233;   // Omega CDM
	cmd.ns=0.97;    // spectral index of linear P(k) 
	cmd.w=-1.0;       // equation of state of dark energy
	cmd.Omnu=0.0;     // Omega_nu  //Massive neutrinos
	
	cmd.Omm=cmd.Omb+cmd.Omc+cmd.Omnu;
	//~ cmd.Omm=0.3;s
	cmd.Omw=1.-cmd.Omm;
	
	cmd.fnamePS="./input/linear_pk_Takahashi_z0.txt";
	
	
	cmd.chatty = 2;
	//~ cmd.tree_level = 2;  // = 1 uses only second order PT
	cmd.prefix = "run1_";	
	cmd.path_Bells = "Bell_outputs";
	cmd.tree_level = 3;  // = 1 second order PT
						 // = 2 : B=pk^
						 // = 3 uses EFT ctr
						 // other Halo model
						 
	
	//~ cmd.zbin = 0.5078; //20/3.;	
	cmd.zbin = 0.5078; //20/3.;	
	cmd.mMax = 5;  // Bm moments upto mMax 
	

	cmd.chiQuadSteps =300;
	cmd.GLpoints =64;
	//~ cmd.chiQuadSteps =20;
	//~ cmd.GLpoints =16;
	
	//~ cmd.Nell = 128;	
	cmd.Nell = 128;	
	cmd.ellmax = 10000.0;
	cmd.ellmin = 0.001;	
	//~ cmd.ellmax = 10000.0;
	//~ cmd.ellmin = 1.0;	
	
	
	cmd.Wg = 0;  // =0 Wg(chi) Dirac Delta at chibin=chi(zBin)
				 // =1 Wg(chi) from file cmd.fWgchi
 
	//~ cmd.fWgchi="./input/Wg_Takahashi_z20548.txt"; 
	cmd.fWgchi="./input/Wg_Takahashi_z05078.txt"; 

	
	if(cmd.chatty>1) printf("OmegaM = %f \n",cmd.Omm);
	cmd.writevectors =1;
	
	
  }






void MainLoop(void)
{

	gv.time=clock();
	background(cmd.zbin);
	if(cmd.chatty>1){
		printf("time evaluating background: %lf s\n", (double)(clock() - gv.time) / CLOCKS_PER_SEC );
		printf("iv.chiMaxInt = %f Mpc/h\n", iv.chiMaxInt);
		printf("gv.chiBin    = %f Mpc/h\n", gv.chiBin);
	 }

	//~ printf("\n chiOfzT= \n {");    
	//~ for	(int i=0;i<gv.Nz;i++) printf("{%f, %f}, ", gv.zT[i], gv.chiOfzT[i]);
	//~ printf("} \n");    
        
	gv.time=clock();
	gv.z=cmd.zbin;
	gv.Dp = Dplusf(gv.z)/gd.Dpz0;
	gv.r_sigma = calcrsigma(gv.Dp, 0.001,8.,100);   
	gv.n_eff= n_eff_func(gv.r_sigma,gv.Dp,0.001,8.,100);  
	BmKspace(cmd.mMax, 0.001, 1.0, 100, cmd.GLpoints, gv.z,gv.Dp,gv.r_sigma,gv.n_eff);    
    printf("time= %lf \n", (double)(clock() - gv.time) / CLOCKS_PER_SEC ); 
     
	allocate_iv();
	gv.time=clock();
	Bmell();
	if(cmd.chatty>1) printf("Integration time= %lf \n", (double)(clock() - gv.time) / CLOCKS_PER_SEC );
	
	allocate_zv();
	get_zetam();    
	

	
	write(); 
	
	//Tests();
	  	
}







void Tests(void)
{	
	//read and write
	//~ FILE *fp;
	//~ int size = iv.Nell*iv.Nell;
	//~ int n_data, ndatamax;
    //~ double Bmell0vector[size];
    //~ double mat[iv.Nell][iv.Nell];
    //~ int i;
    
	//~ fp = fopen("outputs/1_BmellsVector_0.txt", "r");    
     
    //~ ndatamax = size;
    
    //~ n_data=0;
    //~ if (NULL == fp) {
        //~ printf("\n\nfile can't be opened \n\n");
    //~ }
    //~ if(fp!=NULL){   
		//~ while(fscanf(fp, "%lf", &Bmell0vector[n_data])!=EOF){
			//~ n_data++;
			//~ if(n_data>ndatamax) printf("n_data_max should be larger than the number of data lines \n");
		//~ }	  	
	//~ fclose(fp);		
	//~ }
    
 	//~ for(int i=0;i<size; i++){
	//~ printf("%e,  ",Bmell0vector[i]);	
	//~ }   
    
    //~ printf(" \n");
	//~ for(int i=0; i<iv.Nell ; i++){
	//~ for(int j=0; j<iv.Nell ; j++){
			//~ mat[i][j] = Bmell0vector[i*iv.Nell+j];	
	//~ }
	//~ }    
    
    //~ printf(" \n");
	//~ for(int i=0; i<iv.Nell ; i++){
	//~ for(int j=0; j<iv.Nell ; j++){	
			//~ printf("%15e ", mat[i][j]);
			//~ if(j==iv.Nell-1) printf(" \n");
	//~ }
	//~ } 



}














