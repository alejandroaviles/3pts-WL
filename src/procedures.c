#include "functions.h"
#include "procedures.h"

#define m_PI   3.1415926535897932384626433

/////////////////////////////////////////////
/* Initial */
/////////////////////////////////////////////
void Initial()
{
	gd.Dpz0 = Dplusf(0.0);	
	read_inputpk();
	gd.sigma8 =  sigmaRTH(8,0.001,8.,100);	
}


void allocate_iv(void)
{
	iv.Nell   = cmd.Nell;
	iv.ellmax = cmd.ellmax;
	iv.ellmin = cmd.ellmin;
	iv.ellT   = malloc(iv.Nell * sizeof(double *));
	
	for(int i=0; i<iv.Nell ; i++){		
		iv.ellT[i]= exp(log(iv.ellmin) 
		+ i*log(iv.ellmax/iv.ellmin)/(iv.Nell-1.0));	
	}	

	int NumMoments=cmd.mMax+1;
	iv.BmVectors  = malloc(NumMoments * sizeof(double *));
	iv.BmVectorsp = malloc(NumMoments * sizeof(double *));            
	for(int m=0; m<NumMoments; m++) {
		iv.BmVectors [m] = malloc(iv.Nell * iv.Nell * sizeof(double));
		iv.BmVectorsp[m] = malloc(iv.Nell * iv.Nell * sizeof(double));
		for(int ij=0; ij< iv.Nell*iv.Nell; ij++){
			 iv.BmVectors [m][ij] = 0;
			 iv.BmVectorsp[m][ij] = 0;
		}
	} 
	
}


// Needs allocated vectors iv.BmVectors[m]
void Bmell(void)
{	
	double chi, z, Dp, rsigma, neff, qv;
	double chiprev, val, half, deltachi;
	double chimax, chimin;
	int NumMoments=cmd.mMax+1;
	
	double BmvectorsB[NumMoments][iv.Nell*iv.Nell];
	double BmvectorsA[NumMoments][iv.Nell*iv.Nell];
	
	if(cmd.chatty>0){
		printf("\nComputing Bm(ell1,ell2) for symmetric %d x %d array of ell values  \n", 
			iv.Nell, iv.Nell);
		printf("	ellmin = %f, ellmax = %f  \n", 
			iv.ellmin, iv.ellmax);
		printf("	Number of moments = %d\n", cmd.mMax+1);
		printf("	Quadrature chi steps= %d\n", iv.chiQuadSteps);
	}
		
	chimax = iv.chiT_chiint[iv.chiQuadSteps-1];
	chimin = iv.chiT_chiint[0];

	for (int i=0;i<iv.chiQuadSteps;i++){
		chi    = iv.chiT_chiint[i];
		z      = iv.zT_chiint[i];
		Dp     = iv.DpT_chiint[i];
		rsigma = iv.rsigma_chiint[i];
		neff   = iv.neff_chiint[i];
		qv     = iv.q_chiint[i];
		
		deltachi = chi - chiprev;
		
		Bm(chi,z,Dp,rsigma,neff);
		
		for(int m=0; m<NumMoments; m++) {
		for(int ij=0; ij<iv.Nell*iv.Nell ; ij++){
			BmvectorsB[m][ij] = pow(qv,3.)/pow(chi,4.) * iv.BmVectors[m][ij];
		}			
		}	
		
		if(i==0){ 
			for(int m=0; m<NumMoments; m++) {
				for(int ij=0; ij<iv.Nell*iv.Nell ; ij++){
				BmvectorsA[m][ij] = BmvectorsB[m][ij];
				}			
			}			
			deltachi = chi;
		} 
		
		for(int m=0; m<NumMoments; m++) {
			for(int ij=0; ij<iv.Nell*iv.Nell ; ij++){
				iv.BmVectorsp[m][ij] +=  0.5*(BmvectorsA[m][ij]
							  +BmvectorsB[m][ij])* deltachi;						
			}			
		}

		for(int m=0; m<NumMoments; m++) {
		for(int ij=0; ij<iv.Nell*iv.Nell ; ij++){
			BmvectorsA[m][ij] = BmvectorsB[m][ij];
		}			
		}	 
		chiprev = chi;
		
	// TEST	
	//~ FILE *fp;	
	//~ char str[100];
	//~ int ii;
	//~ sprintf(str,"%s/step_%d_Bm0diag.txt",cmd.path_Bells,i);
	//~ fp = fopen (str, "w+");	
		//~ //for(int ij=0; ij<iv.Nell*iv.Nell ; ij++){	
			//~ //fprintf(fp, "%15e \n", BmvectorsA[0][ij]);
		//~ //}	
		//~ for(int kk=0; kk<iv.Nell ; kk++){	
			//~ ii = kk * iv.Nell;
			//~ fprintf(fp, "%15e \n", iv.BmVectors[0][ii]);
		//~ }	
	//~ fclose (fp);	
	//
		
	}

	// TEST	
	FILE *fp;	
	char str[100];
	sprintf(str,"%s/tests/i_z_chi_qv_Dp.txt",cmd.path_Bells);
	fp = fopen (str, "w+");	
		for(int i=0; i<iv.chiQuadSteps ; i++){	
			fprintf(fp, "%d %15e %15e %15e %15e \n", 
			i, iv.zT_chiint[i], iv.chiT_chiint[i], 
			iv.q_chiint[i], iv.DpT_chiint[i]);
		}	
	fclose (fp);	
		
		
}


void Bm(double chi, double z, double Dp, double r_sigma, double n_eff)
{	
	double k1,k2,varphi,w,k3,BT, ell1,ell2;
	double *xGL, *wGL;
	int m=0;
	int NumMoments=cmd.mMax+1;
	double val[NumMoments];
	double EFTctr;
	
	xGL = malloc(cmd.GLpoints * sizeof(double));
	wGL = malloc(cmd.GLpoints * sizeof(double));	
	gaussleg(-m_PI, m_PI, xGL, wGL, cmd.GLpoints);

	EFTctr=-3.0*pow(Dp,2);	

	
	for(int i=0; i<iv.Nell; i++){
	for(int j=i; j<iv.Nell; j++){
		ell1=iv.ellT[i];
		ell2=iv.ellT[j];
		k1 = ell1/chi;
		k2 = ell2/chi;
		for(int m=0;m<NumMoments;m++) val[m] = 0.0;
		
		for(int i=0; i<cmd.GLpoints;i++){
			varphi = xGL[i];
			w   = wGL[i];
			k3  = sqrt( k1*k1 + k2*k2 - 2.*k1*k2 * cos(varphi) );
			if (cmd.tree_level==1){
				BT = Bispec_tree(k1, k2, k3, Dp);
			} else if (cmd.tree_level==2){
				BT = Bispec_P2(k1, k2, k3, Dp);	
			} else if (cmd.tree_level==3){
				BT  = Bispec_tree_EFT(k1, k2, k3, Dp, EFTctr);
			} else {
				BT  = Bispec_Takahashi(k1, k2, k3, z, Dp, r_sigma, n_eff);
			}
			for(int m=0;m<NumMoments;m++) val[m] =  val[m] + w*BT*cos(m*varphi);			
		}

		for(int m=0; m<NumMoments; m++){
			iv.BmVectors[m][i*iv.Nell + j] = val[m]/(2*m_PI);
			if(j!=i) iv.BmVectors[m][j*iv.Nell + i] = val[m]/(2*m_PI);
		}	
	}
	}				
				
}


void BmKspace(int Maxm, double kmin, double kmax, int Nk, int GLpoints, double z, double Dp, double r_sigma, double n_eff)
{	
	if(cmd.chatty>0) printf("\nComputing Bm(k1,k2) for symmetric array of %d x %d k1,k2 values  \n", Nk, Nk);	
	
	double k1,k2,varphi,w,k3,BT;
	double *kT,*xGL, *wGL;
	double mat[Maxm+1][Nk][Nk];
	double val[Maxm+1];
	int m=0;

	kT  = malloc(Nk       * sizeof(double));
	xGL = malloc(GLpoints * sizeof(double));
	wGL = malloc(GLpoints * sizeof(double));
	
	gaussleg(-m_PI, m_PI, xGL, wGL, GLpoints);
    
	if(cmd.chatty==2) 	printf("Maxm=%d, kmin=%f ,kmax=%f, Nk=%d, GLpoints=%d \n", 
	        Maxm,    kmin,    kmax,    Nk,    GLpoints);
	if(cmd.chatty==2) printf("z=%f, Dp=%f, rsigma=%f, neff=%f \n", z,Dp,r_sigma,n_eff);         
	
	for(int i=0; i<Nk ; i++){		
		kT[i]= exp(log(kmin) + i*log(kmax/kmin)/(Nk-1.0));	
	}

	gv.time=clock();

	for(int i=0; i<Nk ; i++){
	for(int j=i; j<Nk ; j++){
		k1=kT[i];
		k2=kT[j];
		for(int m=0;m<=Maxm;m++) val[m] = 0.0;
		
		for(int i=0; i<GLpoints;i++){
			varphi = xGL[i];
			w   = wGL[i];
			k3  = sqrt( k1*k1 + k2*k2 - 2.*k1*k2 * cos(varphi) );
			BT  = Bispec_Takahashi(k1, k2, k3, z, Dp, r_sigma, n_eff);
			for(int m=0;m<=Maxm;m++) val[m] =  val[m] + w*BT*cos(m*varphi);			
		}

		for(int m=0; m<=Maxm; m++){
			mat[m][i][j] = val[m];
			if(j!=i) mat[m][j][i]=val[m];
		}
		
	}
	}	
	
	if(cmd.chatty==2) printf("Bnk time = %f15 \n", (double)(clock()-gv.time) / CLOCKS_PER_SEC  );
	if(cmd.chatty==2) printf("\n");
	
	// Write
	char int_str[20];
	char str[100];
	for (int m=0; m<Maxm+1; m++){
		FILE *fp;
		sprintf(int_str, "%d", m);
		sprintf(str,"%s/%sBnk_%s.txt",cmd.path_Bells,cmd.prefix,int_str);		
		fp = fopen (str, "w+");
		
		for(int i=0; i<Nk ; i++){
		for(int j=0; j<Nk ; j++){	
			fprintf(fp, "%15e   ", mat[m][i][j]);
			if(j==Nk-1 && i!=Nk-1) fprintf(fp, " \n");
		}
		}
		
		fclose (fp);
	}	

	FILE *fp;
	sprintf(str,"%s/%skArray.txt",cmd.path_Bells,cmd.prefix);	
	fp = fopen (str, "w+");		
	for(int i=0; i<Nk ; i++){		
		fprintf(fp,"%15e\n", kT[i]);		
	}
	fclose (fp);	

}



/////////////////////////////////////////////
/* Routines */
/////////////////////////////////////////////
#define EPS 3.0e-11
void gaussleg(double x1, double x2, double xGL[], double wGL[], int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;
    
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++){
		z=cos(m_PI*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		xGL[i-1]=xm-xl*z;
		xGL[n+1-i-1]=xm+xl*z;
		wGL[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		wGL[n+1-i-1]=wGL[i-1];
	}
}
#undef EPS


// Read and extrapolate input pk
void read_inputpk() //Extrapolation not yet implemented
{
	FILE *fp;	
	gd.n_data=0;
	fp=fopen(cmd.fnamePS,"r");   // linear P(k) table
	
    if (NULL == fp) {
        printf("\n\nlinear power spectrum can't be opened \n\n");
    }
	
	if(fp!=NULL){   // input: k[h/Mpc]   P(k)[(Mpc/h)^3]
		while(fscanf(fp, "%lf %lf", &gd.k_data[gd.n_data], &gd.pkz0_data[gd.n_data])!=EOF){
			gd.n_data++;
			if(gd.n_data>n_data_max) printf("n_data_max should be larger than the number of data lines \n");
		}	  	
	fclose(fp);		
	}		
}



double interpolation1(double x, double xT[], double yT[], int n_data)   // interpolation order 1
{
  int j,j1,j2,jm;
  double f;

  if(x<xT[0]) return 0.;
  if(x>xT[n_data-1]) return 0.;
  
  j1=0, j2=n_data-1, jm=(j1+j2)/2;
  for(;;){
    if(x>xT[jm]) j1=jm;
    else j2=jm;
    jm=(j1+j2)/2;

    if(j2-j1==1) break;
  }
  j=j1;

  f=(yT[j+1]-yT[j])/(xT[j+1]-xT[j]) * (x -xT[j])+ yT[j];
  
  return f;  
}


double interpolationlog(double x, double xT[], double yT[], int n_data)   // interpolation in log space
{
  int j,j1,j2,jm; 
  double f;

  if(x<xT[0]) return 0.;
  if(x>xT[n_data-1]) return 0.;
  
  j1=0, j2=n_data-1, jm=(j1+j2)/2;
  for(;;){
    if(x>xT[jm]) j1=jm;
    else j2=jm;
    jm=(j1+j2)/2;

    if(j2-j1==1) break;
  }
  j=j1;

  f=(log10(yT[j+1])-log10(yT[j]))/(log10(xT[j+1])
      -log10(xT[j]))*(log10(x)-log10(xT[j]))+log10(yT[j]);
  
  return pow(10.,f);  
}



/////////////////////////////////////////////
/* Write */
/////////////////////////////////////////////


void write(void)
{
	double elli, ellj;	
	char int_str[20];
	char str[100];	
	
	for (int m=0; m<cmd.mMax+1; m++){
		FILE *fp;
		sprintf(int_str, "%d", m);
		sprintf(str,"%s/%sBmells_%s.txt",cmd.path_Bells,cmd.prefix,int_str);
		sprintf(int_str, "%d", m);
		fp = fopen (str, "w+");
		
		for(int i=0; i<iv.Nell ; i++){
			for(int j=0; j<iv.Nell ; j++){	
				fprintf(fp, "%15e ", iv.BmVectorsp[m][i*iv.Nell + j]);
				if(j==iv.Nell-1 && i!=iv.Nell-1) fprintf(fp, " \n");
			}
		}		
		fclose (fp);
	}	
		
	
	FILE *fp2;	
	sprintf(str,"%s/%sellArray.txt",cmd.path_Bells,cmd.prefix);
	fp2 = fopen (str, "w+");		
	for(int i=0; i<iv.Nell ; i++){		
		fprintf(fp2,"%15e\n", iv.ellT[i]);		
	}
	fclose (fp2);
		
	FILE *fp3;
	sprintf(str,"%s/%sinfo.txt",cmd.path_Bells,cmd.prefix);
	fp3 = fopen (str, "w+");		
	fprintf(fp3, "Cosmological Parameters: \n");
	fprintf(fp3, "	Omega_m = %f\n",cmd.Omm);
	fprintf(fp3, "	     ns = %f\n",cmd.ns);
	fprintf(fp3, "\n");
	fprintf(fp3, "Computing Bm(ell1,ell2) for %d x %d array  \n", 
			iv.Nell, iv.Nell);
	fprintf(fp3, "	ellmin = %f, ellmax = %f  \n", 
			iv.ellmin, iv.ellmax);
	fprintf(fp3, "	Number of moments = %d\n", cmd.mMax+1);
	fprintf(fp3, "	Quadrature chi steps= %d\n", iv.chiQuadSteps);
	fclose (fp3);

	if (cmd.writevectors ==1){
	for (int m=0; m<cmd.mMax+1; m++){
		FILE *fp;
		sprintf(int_str, "%d", m);
		sprintf(str,"%s/%sBmellsVector_%s.txt",cmd.path_Bells,cmd.prefix,int_str);
		sprintf(int_str, "%d", m);
		fp = fopen (str, "w+");
		
		for(int ij=0; ij<iv.Nell*iv.Nell ; ij++){	
			fprintf(fp, "%e \n", iv.BmVectorsp[m][ij]);
		}		
		fclose (fp);
	}
	}

}




/////////////////////////////////////////////
/* free variables */
/////////////////////////////////////////////

//~ void free_variables(void)
//~ {	
		//~ for (int m=0; m<cmd.mMax+1; m++){
			//~ free(iv.BmVectorsp[m]);
			//~ free(iv.BmVectors [m]);
		//~ }	
		//~ free(iv.BmVectorsp);
		//~ free(iv.BmVectors );	
//~ }




















