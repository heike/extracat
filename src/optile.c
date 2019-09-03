#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 



#define MAT_ELT(x, i, j, nrow) x[(long)(i)+(long)(j)*(long)(nrow)]


int mymax(int a, int b){
	if(a < b){
		return b;
	}else{
		return a;
	}
}
int mymin(int a, int b){
	if(a < b){
		return a;
	}else{
		return b;
	}
}

float mymin2(float a, float b){
	if(a < b){
		return a;
	}else{
		return b;
	}
}

float mymax2(float a, float b){
	if(a < b){
		return b;
	}else{
		return a;
	}
}

// ----------------------------------------------------------------------------------

static void quickSort (float *a, int *order, int lo, int hi)
{
	//  lo is the lower index, hi is the upper index
	//  of the region of array a that is to be sorted
    int i=lo, j=hi, ho, hi2,lo2;
	float atmp;
	
    float x=a[(lo+hi)/2];
	//	Rprintf("from %d to %d with pivot a[%d] = %f\n\n",i,j,(lo+hi)/2,x);
    //  partition
    do
    {    
        while (a[i]<x) i++; 
        while (a[j]>x) j--;
        if (i<=j)
        {
            atmp=a[i]; a[i]=a[j]; a[j]=atmp;
			ho = order[i]; order[i]=order[j]; order[j]=ho;
            i++; j--;
			//Rprintf("ex %d -- %d\n\n",i,j);
        }
    } while (i<=j);
	hi2 = j;
	lo2 = i;
    //  recursion
    if (lo<j) quickSort(a, order, lo, hi2);
    if (i<hi) quickSort(a, order, lo2, hi);
	
}




// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //





int getindex(int dims[], int ind[], int nd){
	
	int i,s;
	//int nd2 = sizeof(dims)/2;
	for (i=0; i<nd; i++) {
		if (ind[i] > dims[i]-1) {
			Rprintf("invalid index!\n");
			//for (int i2=0; i2<nd; i2++) {
			//	Rprintf("%d vs. %d,\t",ind[i2],dims[i2]);
			//}
			return -1;
		}
	}
	int cpx[nd];
	cpx[0] = dims[0];
	
	for (i = 1; i < nd; i++) {
		cpx[i] = cpx[i-1]*dims[i];
	}
	
	
	int k = ind[0]; 
	
	for (s = 1; s < nd; s++) {
		//Rprintf(" %d + %d * %d\n",k,ind[s],cpx[s-1]);
		k = k + ind[s]*cpx[s-1];
	}
	return k;
} 




// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //



void diagprod(int *mv1, int *mv0, float *cv, int *IM, int dims[], int index[], int step, int nd){
	int i,j, s;
	
	//int nd = sizeof(dims)/2;
	int val1, val2;
	
	int (*IX)[nd] = (int (*)[nd])IM;
	float tmp, tmp2;
	//Rprintf("\n welcome to diagprod\n");
	
	//Rprintf("step = %d\n",step);
	//Rprintf("nd = %d\n",nd);
	int indS[nd];
	for (i = 0; i< nd; i++) {
		indS[i] = index[i];
	}
	
	
	if( step == nd-1 ){ // the recursion in in the last dimension
		int ind2[nd];
		//Rprintf("\n IX in diagprod:\n");
		//for (int r = 0; r < nd; r++) {
		//	for (int t = 0; t < dims[r]; t++) {
		//		Rprintf("%d\t",IX[t][r]);
		//	}
		//	Rprintf("\n");
		//}
		for (j = 0; j < nd-1; j++) {
			//Rprintf("ind %d is %d\n",j,indS[j]);
			ind2[j] = IX[indS[j]-1][j];
			indS[j] = IX[indS[j]][j];
			//Rprintf("and changes to %d\t",indS[j]);
		}
		for (s = 1; s < dims[nd-1]; s++) {
			indS[nd-1] = IX[s][nd-1];
			ind2[nd-1] = IX[s-1][nd-1];
			//Rprintf("final ind:\n");
			//for (i=0; i<nd; i++) {
			//	Rprintf("%d\t",indS[i]);
			//}
			//Rprintf("\n");
			//for (i=0; i<nd; i++) {
			//	Rprintf("%d\t",ind2[i]);
			//}
			//Rprintf("\n ind1 = %d",getindex(dims,indS,nd));
			//Rprintf("\n ind2 = %d",getindex(dims,ind2,nd));
			//Rprintf("\n val1 = %d",mv0[ getindex(dims,indS,nd) ]);
			//Rprintf("\n val2 = %d",mv1[ getindex(dims,ind2,nd) ]);
			val1 = mv0[ getindex(dims,indS,nd) ];
			val2 = (int) mv1[ getindex(dims,ind2,nd) ];
			tmp = (float) val1*val2;
			tmp2 = cv[0];
			//Rprintf("\n tmp = %f",tmp);
			cv[0] = tmp + tmp2;//((float*) val1) * ((float*) val2);
			//Rprintf("\n cv = %f",cv[0]);
		}
	}else {
		
		for (s = 1; s < dims[step]; s++) {
			indS[step] = s;
			//Rprintf("temp ind:\n");
			//for (i=0; i<nd; i++) {
			//	Rprintf("%d\t",indS[i]);
			//}
			//Rprintf("\n");
			diagprod(mv1,mv0, cv, IM, dims, indS, step+1,  nd);
		}
	}
	
}





// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //


float mvclasscrit2(int *m0, int *IM, int dims[], int nd){
	
	
	int i, j1,j2, s1, s2, k, rind ,lind;
	//int nd = sizeof(dims)/2;
	//Rprintf("welcome to classcrit\n");
	//Rprintf("dim = %d\n",nd);
	
	int (*IX)[nd] = (int (*)[nd])IM;
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	//for (i=0; i<nd; i++) {
	//	Rprintf("%d\t",dims[i]);
	//}
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dims[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dims[nd-i];//CHECK
	}
	//for (i = 0; i < nd; i++) {
	//	Rprintf("%d\t",cp[i]);
	//}
	//for (i = 0; i < nd; i++) {
	//	Rprintf("%d\t",cpr[i]);
	//}
	// m1 array for the cumsums
	int ml = cp[nd-1]*dims[nd-1];
	int m1[ ml ];
	//int m2[ ml ];
	for (i = 0; i < ml; i++) {
		m1[i] = m0[i];
		//m2[i] = m0[i];
		//Rprintf("%d\t",m1[i]);
	}
	//Rprintf("m1 raw:\n");
	//for (i=0; i<ml; i++) {
	//	Rprintf("%d\t",m1[i]);
	//}
	//go through m1 and m0 with stepsize s acc. to dimension k
	// dims is in reverse order
	for (k = 0; k < nd; k++) {
		s1 = cp[k];
		s2 = cpr[k];
		//Rprintf("s1 = %d\n",s1);
		//Rprintf("s2 = %d\n",s2);
		//if (dims[k] > 2) {
		for (i = 1; i < dims[k]; i++) {
			for (j1 = 0; j1 < s2; j1++) {
				for (j2 = 0; j2 < s1; j2++) {
					rind = IX[i][k]*s1 + j1*s1*dims[k] + j2;
					lind = IX[i-1][k]*s1 + j1*s1*dims[k] + j2;
					//Rprintf("%d\t",rind );
					//Rprintf("%d\n",lind);
					m1[ rind ] += m1[ lind ];
				}
			}
		}
		//m1[  j1*s1*dims[k] + j2 - 1] = m0[ j1*s1*dims[k] + j2 - 1];
		//} 
		
	}
	
	// multiplication of m0 with m1[ ind -1 ] via diagprod
	float cv = 0.0;
	int *m1p = &m1[0];
	float *cvp = &cv;
	int init[nd];
	for (i=0; i<nd; i++) {
		init[i] = 0;
	}
	diagprod(m1p,m0,cvp,IM,dims,init,0, nd);
	//Rprintf("crit = %f\n",cv);
	return cv;
}


// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //


SEXP simmich(SEXP P, SEXP Q, SEXP n, SEXP N){
	int i,j;
	int m = INTEGER(n)[0];
	int M = INTEGER(N)[0];
	//Rprintf("%d, %d\n\n",M,m);
	double *Pv = calloc(M,sizeof(double));
	double *Qv = calloc(M,sizeof(double));
	
	double sp = 0.0;
	for (i = 0; i < m; i++) {
		Pv[i] = REAL(P)[i];
		Qv[i] = REAL(Q)[i];
	}
	//Rprintf("%f\n",M_PI);
	//Rprintf("%f\n\n\n",log(M_PI));
	for (i = m; i < M; i++) {
		Pv[i] = 1/(M_PI*i);
		Qv[i] = Pv[i];
		for (j = 0; j < i; j++) {
			Qv[i] -= Qv[j]*Pv[i-j-1];
		}
		//Rprintf("%f\n",Qv[i]);
	}
	//Rprintf("%f, %f",Pv[M-1],Qv[M-1]);
	double b = 0.0;
	for (i = 0; i < M; i++) {
		b += Qv[i] * log(i+1);
		sp += Qv[i];
	}
	//Rprintf("%f\n\n\n",b);
	b -= M_PI*log(log(M));
	
	SEXP ret = allocVector(REALSXP,2);
	REAL(ret)[0] = b ;
	REAL(ret)[1] = sp ;
	free(Pv);
	free(Qv);
		
	return ret;

}


SEXP hammdist(SEXP dset){
	
	int n;
	int m;
	int i, i2;
	int j,k,d,d2;
	
	int dist;
	
	SEXP sDim = getAttrib(dset, R_DimSymbol);
	Rprintf("test");
	
	if (isInteger(sDim) && (LENGTH(sDim) == 2)) {
		int *dim = INTEGER(sDim);
		Rprintf("%d x %d matrix\n", dim[0], dim[1]);
		n = dim[0];
		m = dim[1];
	} else error("invalid dimensions");
	
	SEXP dm = allocVector(INTSXP,n*(n-1)/2);
	k = 0;
	for( i=0; i < n-1; i++ ){
		for( i2 = i+1; i2 < n; i2++ ){
			
			dist = 0;
			
			
			for( j=0; j < m; j++ ){
				d = n*j + i;
				d2 = n*j + i2;
				dist = (INTEGER(dset)[d] == INTEGER(dset)[d2]) ? dist : dist+1;
			}
			INTEGER(dm)[k] = dist;
			k = k+1;
			
		}
	}
	
	
	return dm;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP dcorD(SEXP x, SEXP y, SEXP freq){
	
	
	int n = LENGTH(freq);
	int s = n*(n-1)/2;
	int i, j,k;
	
	
	double *DMY = calloc(s,sizeof(double));
	double *DMX = calloc(s,sizeof(double));
	double *F = calloc(s,sizeof(double));
	
	double Edx[n];
	double Edy[n];
	for (i=0; i<n; i++) {
		Edx[i] = 0;
		Edy[i] = 0;
	}
    double S1 = 0;
    double S2 = 0;
    double S3 = 0;
    
    double S2a = 0;
    double S2b = 0;
    
    double  S1X = 0;
	double  S1Y = 0;
    double  S2X = 0;
	double  S2Y = 0;
	double  S3X = 0;
	double  S3Y = 0;
	
	k = 0;
	for (i=0; i<(n-1); i++) {
		for (j = i+1; j < n; j++) {
			DMX[k] = REAL(x)[i+j*n];
			DMY[k] = REAL(y)[i+j*n];
		
			F[k] = REAL(freq)[i] * REAL(freq)[j];
			S1 += DMX[k]*DMY[k]*F[k];
			S1X += DMX[k]*DMX[k]*F[k];
			S1Y += DMY[k]*DMY[k]*F[k];
			
			Edx[i] += DMX[k]*REAL(freq)[j];
			
			Edy[j] += DMY[k]*REAL(freq)[i];
			Edx[j] += DMX[k]*REAL(freq)[i];
			Edy[i] += DMY[k]*REAL(freq)[j];
			k++;
		}
	}
	
	
	for (i = 0; i < n; i++) {
		S3 += Edx[i] * Edy[i] * REAL(freq)[i];
		S2a += Edy[i] * REAL(freq)[i]; 
		S2b += Edx[i] * REAL(freq)[i];
		S3X += Edx[i] * Edx[i] * REAL(freq)[i];
		S3Y += Edy[i] * Edy[i] * REAL(freq)[i];
	}
	S1 = 2*S1;
	S1Y = 2*S1Y;
	S1X = 2*S1X;
	S2 = S2a*S2b;
	S2X = S2b*S2b;
	S2Y = S2a*S2a;
	
	SEXP ret = allocVector(REALSXP,1);
	REAL(ret)[0] = pow( (S1+S2-2*S3)/pow( (S1X+S2X-2*S3X)*(S1Y+S2Y-2*S3Y) ,0.5) , 0.5) ;
	
	free(DMY);
	free(DMX);
	free(F);
	
	
	return ret;
}





// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //




// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP dcorR(SEXP x, SEXP y, SEXP freq, SEXP e){
	
	
	int n = LENGTH(y);
	int s = n*(n-1)/2;
	
	int i, j,k;
	
		
	double *DMY = calloc(s,sizeof(double));
	double *DMX = calloc(s,sizeof(double));
	double *F = calloc(s,sizeof(double));
	
	double Edx[n];
	double Edy[n];
	for (i=0; i<n; i++) {
		Edx[i] = 0;
		Edy[i] = 0;
	}
	double S1 = 0;
    double S2 = 0;
    double S3 = 0;
    
    double S2a = 0;
    double S2b = 0;
    
    double  S1X = 0;
	double  S1Y = 0;
    double  S2X = 0;
	double  S2Y = 0;
	double  S3X = 0;
	double  S3Y = 0;
	
	k = 0;
		for (i=0; i<(n-1); i++) {
			for (j = i+1; j < n; j++) {
				DMX[k] = pow(fabs(REAL(x)[i] - REAL(x)[j]),REAL(e)[0]);
				DMY[k] = pow(fabs(REAL(y)[i] - REAL(y)[j]),REAL(e)[0]);
              
				F[k] = REAL(freq)[i] * REAL(freq)[j];
                //Rprintf("\t%f, ",F[k]);
				S1 += DMX[k]*DMY[k]*F[k];
				S1X += DMX[k]*DMX[k]*F[k];
				S1Y += DMY[k]*DMY[k]*F[k];
		
				Edx[i] += DMX[k]*REAL(freq)[j];
				//Rprintf("\n%f, ",Edx[i]);
				Edy[j] += DMY[k]*REAL(freq)[i];
				Edx[j] += DMX[k]*REAL(freq)[i];
				Edy[i] += DMY[k]*REAL(freq)[j];
				k++;
			}
           // Rprintf("\n");
		}
	
	//Rprintf("\nEdx = \n");
	//for (i = 0; i < n; i++) {
	//	Rprintf("%f, ",Edx[i]);
	//}
	//Rprintf("\nEdy = \n");
	//for (i = 0; i < n; i++) {
	//	Rprintf("%f, ",Edy[i]);
	//}
	
		for (i = 0; i < n; i++) {
			S3 += Edx[i] * Edy[i] * REAL(freq)[i];
			S2a += Edy[i] * REAL(freq)[i]; 
			S2b += Edx[i] * REAL(freq)[i];
			S3X += Edx[i] * Edx[i] * REAL(freq)[i];
			S3Y += Edy[i] * Edy[i] * REAL(freq)[i];
		}
		S1 = 2*S1;
		S1Y = 2*S1Y;
		S1X = 2*S1X;
		S2 = S2a*S2b;
		S2X = S2b*S2b;
		S2Y = S2a*S2a;
		
		
	//Rprintf("S1 = %f\nS2 = %f\nS3 = %f\n",S1,S2,S3);
	//Rprintf("S1X = %f\nS2X = %f\nS3X = %f\n",S1X,S2X,S3X);
	//Rprintf("S1Y = %f\nS2Y = %f\nS3Y = %f\n",S1Y,S2Y,S3Y);
	
	SEXP ret = allocVector(REALSXP,1);
	REAL(ret)[0] = pow( (S1+S2-2*S3)/pow( (S1X+S2X-2*S3X)*(S1Y+S2Y-2*S3Y) ,0.5) , 0.5) ;
	//Rprintf("%f, %f, %f, %f, %f, %f, %f, %f, %f",S1,S2,S3,S1X,S2X,S3X,S1Y,S2Y,S3Y);
	
	free(DMY);
	free(DMX);
	free(F);
	
	//Rprintf("S1 = %f, S2 = %f, S3 = %f,S1X = %f,S2X = %f,S3X = %f\n\n",S1,S2,S3,S1X,S2X,S3X);

	return ret;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP dcorM(SEXP x, SEXP dim, SEXP freq, SEXP e){
	
	int m = INTEGER(dim)[0];
	int n = LENGTH(x)/m;
	int s = n*(n-1)/2;
	int sm = s*m;
	int nm = n*m;
	int i, j;
	int k = 0;
	int vi, vj;
	
	//Rprintf("m = %d, n = %d , s = %d, e = %f\n\n---",n,m,s,REAL(e)[0]);
	
	double *Edx = calloc(nm,sizeof(double));
	
	// distances
	double *DMX = calloc(sm,sizeof(double));
	
	// frequencies/weights, 1 per pair
	double *F = calloc(s,sizeof(double));


    
	//for (i=0; i<n; i++) {
	//	for (vi = 0 ; vi < m; vi++) {
	//		MAT_ELT(Edx,i,vi,nm) = 0;
	//	}
	//}
	
	
	double *S1 = calloc(m*m,sizeof(double));
	double *S2 = calloc(m*m,sizeof(double));
	double *S3 = calloc(m*m,sizeof(double));
	
	double *S1X = calloc(m,sizeof(double));
	double *S2X = calloc(m,sizeof(double));
	double *S3X = calloc(m,sizeof(double));
	
	double *S2ab = calloc(m,sizeof(double));
	
	
		
	
	for (i=0; i<(n-1); i++) {
				for (j = i+1; j < n; j++) {//k
		
					
					F[k] = REAL(freq)[i] * REAL(freq)[j];
					
					for (vi = 0; vi < m; vi++) {
                        // WRITING to lower.tri
                        
					
						MAT_ELT(DMX,k,vi,s) = pow(fabs(REAL(x)[i + vi*n] - REAL(x)[j+vi*n]),REAL(e)[0]);
						
						S1X[vi] += MAT_ELT(DMX,k,vi,s)*MAT_ELT(DMX,k,vi,s)*F[k];
						
                        MAT_ELT(Edx,i,vi,n) += MAT_ELT(DMX,k,vi,s)*REAL(freq)[j];
						MAT_ELT(Edx,j,vi,n) += MAT_ELT(DMX,k,vi,s)*REAL(freq)[i];
                        
                       if(vi > 0){
                            for (vj = 0; vj < vi; vj++) {
                                S1[vi + vj*m] += MAT_ELT(DMX,k,vi,s)*MAT_ELT(DMX,k,vj,s)*F[k];
                            }
                        }
                        
					}
					
                    k++;
                }
                
			}
    
    for (vi = 0; vi < m; vi++) {
        for (i = 0; i < n; i++) {
            S2ab[vi] += MAT_ELT(Edx,i,vi,n) * REAL(freq)[i];
            S3X[vi] += MAT_ELT(Edx,i,vi,n)*MAT_ELT(Edx,i,vi,n) * REAL(freq)[i];
        }
        S1X[vi] = 2*S1X[vi];
        S2X[vi] = S2ab[vi]*S2ab[vi];
        
    }
    
    SEXP ret = allocVector(REALSXP,m*(m-1)/2);
	
    k = 0;
    for (vj = 0; vj < m-1; vj++) {
        for (vi = vj+1; vi < m; vi++) {
            for (i = 0; i < n; i++) {
				S3[vi+vj*m] += MAT_ELT(Edx,i,vi,n)*MAT_ELT(Edx,i,vj,n) * REAL(freq)[i];
            }
           
            S2[vi+vj*m] = S2ab[vi]*S2ab[vj];
            S1[vi+vj*m] = 2*S1[vi+vj*m];
            
			REAL(ret)[k] = pow( (S1[vi+vj*m]+S2[vi+vj*m]-2*S3[vi+vj*m])/pow( (S1X[vi]+S2X[vi]-2*S3X[vi])*(S1X[vj]+S2X[vj]-2*S3X[vj]) ,0.5) , 0.5) ;
            k++;
        }
    }
    
    
			//REAL(ret)[0] = pow( (S1+S2-2*S3)/pow( (S1X+S2X-2*S3X)*(S1Y+S2Y-2*S3Y) ,0.5) , 0.5) ;
		
    free(F);
    free(DMX);
	
	free(Edx);
	free(S1);
    free(S2);
    free(S3);
    
    free(S1X);
    free(S2X);
    free(S3X);
    
    free(S2ab);
	
	return ret;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //




int bincoef(int n, int k){
	int s = 1;

	for(int i = (k+1); i < n+1; i++){
		s = s*i/(i-k);
	}
	
	return s;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP MEsimple(SEXP M, SEXP dims){
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int i,j;
	
	 
	
	
	
	float ME = 0.0;
	
	
	
	//set the first row and the last column
	
	for( i = 0; i < n-1; i++ ){ 
		for (j=0; j<m; j++) {
			
			ME += REAL(M)[n*j + i] * REAL(M)[n*j + i+1]; 
		}
	}
	
	for( i = 0; i < n; i++ ){ 
		for (j=0; j<m-1; j++) {
			
				ME += REAL(M)[n*j + i] * REAL(M)[n*(j+1) + i]; 
		}
	}	
	
	
	
	SEXP out = allocVector(REALSXP,1);
	REAL(out)[0] = ME;
	return out;
	
	
}


// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//




// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP simplecrit(SEXP M, SEXP dims){
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int i,j;
	
	int MX[n][m];
	long MY[n][m];
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			MX[i][j] = INTEGER(M)[n*j + i]; 
		}
	}
	
	
	
	float loss = 0.0;
	
	
	
	//set the first row and the last column
	
	for( i = 0; i < n; i++ ){ 
		MY[i][m-1] = MX[i][m-1];
	}
	
	for( j = 1; j < m; j++ ){ 
		MY[0][j] = MX[0][j];
	}
	
	for( j= m - 2; j > 0 ; j-- ){
		for( i=0; i < n-1; i++ ){
			MY[i][j] = MX[i][j] + MY[i][j+1];
		}
	}
	for( i=1; i < n-1; i++ ){
		for( j = 1 ; j < m  ; j++ ){
			MY[i][j] = MX[i][j] + MY[i-1][j];
			// add all rows but the first
			loss = loss + MX[i+1][j-1]*MY[i][j];
		}
	}
	
	// add the first row terms
	for( j = 1 ; j < m  ; j++ ){
		loss = loss + MX[1][j-1]*MY[0][j];
	}
	
	
	
	
	
	SEXP out = allocVector(REALSXP,1);
	REAL(out)[0] = loss;
	return out;
	
	
}


// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//



// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//


SEXP barysort(SEXP M, SEXP dims, SEXP pv , SEXP vs){
	
	int i,j;//i2, j2, ki, kj, rd, cd
	
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	int nd = 2;
	
	
	//int MX[n][m];
	int nm = mymax(n,m);
	int IX[nm][2];
	int TX[nm][2];
	
	float rw[n];
	float cw[m];
	
	
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	int *mv = &INTEGER(M)[0];
    //int *IXp = &IX[0][0];
	int *TXp = &TX[0][0];
	
	int tix[n];
	int tjx[m];
	
	
	
    //float* MX = calloc(n*m,sizeof(float));
    // one unnecessary copy here...
	
	float v[nm];
	//int no[nm];
	//int no2[nm];
	
	float* iv = &v[0];
	int* ptix = &tix[0];
	int* ptjx = &tjx[0];
	
	for (i = 0; i < n; i++) {
		rw[i] = 0;
	}
	
	for (j = 0; j < m; j++) {
		//MAT_ELT(IX,j,1,n) = j;
		//MAT_ELT(TX,j,1,n) = j;
		IX[j][1] = j;
		TX[j][1] = j;
		tjx[j] = j;
		cw[j] = 0;
		for (i = 0; i < n; i++) {
			//MAT_ELT(IX,i,2,n) = i;
			//MAT_ELT(TX,i,2,n) = i;
			IX[i][0] = i;
			TX[i][0] = i;
			tix[i] = i;
			//MAT_ELT(MX,i,j,n) = REAL(M)[i+j*n];//INTEGER
			cw[j] += mv[i+j*n];   //MAT_ELT(MX,i,j,n);
			rw[i] += mv[i+j*n]; //MAT_ELT(MX,i,j,n);
		}
	}
	int opt = 0;
	float crt0, crt1;
	float crt = 0;
	int row = 0;
	
	while(opt == 0){
		opt = 1;
		for (i=0; i < n; i++) {
			v[i] = 0;
			for (j = 0; j < m; j++) {
				//v[i] += MAT_ELT(MX,   MAT_ELT(IX,i,1,n), MAT_ELT(IX,j,2,n),n)*j;
				v[i] += mv[tix[i]+tjx[j]*n]*j;   //MAT_ELT(MX,   tix[i], tjx[j],n)*j;
			}
			v[i] = v[i]/rw[tix[i]];
		}

		
		
					
		
		// get new row order
		quickSort(iv,ptix,0,n-1);
		

		
		for (j=0; j < m; j++) {
			v[j] = 0;
			for (i = 0; i < n; i++) {
				v[j] += mv[tix[i]+tjx[j]*n]*i;//MAT_ELT(MX,   tix[i], tjx[j],n)*i;
			}
			v[j] = v[j]/cw[tjx[j]];
		}
			
		// get new column order given the new row order
		quickSort(iv,ptjx,0,m-1);
		
				
		for (i = 0; i < n; i++) {
			TX[i][0] = tix[i];
		}
		// compute BCC after row order changes
		crt0 = mvclasscrit2( mv, TXp, dimv, nd);
		
		for (j = 0; j < m; j++) {
			TX[j][1] = tjx[j];
		}
		// compute BCC after row AND column order changes
		crt1 = mvclasscrit2( mv, TXp, dimv, nd);
		
		
		
		if (crt0 > crt) {
			opt = 0;
			crt = crt0;
			for (i=0; i<n; i++) {
				IX[i][0] = tix[i];
			}
			row = 1;
		}
		if (crt1 > crt) {
			opt = 0;
			crt = crt1;
			for (j=0; j<m; j++) {
				IX[j][1] = tjx[j];
			}
			row = 0;
		}
		
		//if (crt1 > crt | crt0 > crt) {
			//better than before... keep it and go on.
		//	opt = 0;
		//	crt = crt1;
		//	for (j=0; j<m; j++) {
				//MAT_ELT(IX,j,1,n) = MAT_ELT(TX,j,1,n);
				//MAT_ELT(IX,j,2,n) = MAT_ELT(TX,j,2,n);
		//		IX[j][1] = tjx[j];
		//	}
		//	for (i=0; i<n; i++) {
				//MAT_ELT(IX,j,1,n) = MAT_ELT(TX,j,1,n);
				//MAT_ELT(IX,j,2,n) = MAT_ELT(TX,j,2,n);
				//IX[i][0] = tix[i];
		//	}
		//}else {
			// do nothing and break;
		//}
		if( (crt1 < crt) & (crt0 < crt) ){
			// finish
			if(row == 1){
				for (i = 0; i < n; i++) {
					tix[i] = IX[i][0];
				}
				for (j=0; j < m; j++) {
					v[j] = 0;
					for (i = 0; i < n; i++) {
						v[j] += mv[tix[i]+tjx[j]*n]*i;//MAT_ELT(MX,   tix[i], tjx[j],n)*i;
					}
					v[j] = v[j]/cw[tjx[j]];
				}
				
				// get new column order given the new row order
				quickSort(iv,ptjx,0,m-1);
				for (j=0; j<m; j++) {
					IX[j][1] = tjx[j];
				}
			}else{
				// get best column order
				for (j = 0; j < m; j++) {
					tjx[j] = IX[j][1];
				}
				for (i=0; i < n; i++) {
					v[i] = 0;
					for (j = 0; j < m; j++) {
						//v[i] += MAT_ELT(MX,   MAT_ELT(IX,i,1,n), MAT_ELT(IX,j,2,n),n)*j;
						v[i] += mv[tix[i]+tjx[j]*n]*j;   //MAT_ELT(MX,   tix[i], tjx[j],n)*j;
					}
					v[i] = v[i]/rw[tix[i]];
				}
				// reorder rows
				quickSort(iv,ptjx,0,m-1);
				for (i=0; i<n; i++) {
					IX[i][0] = tix[i];
				}
			}
		}
		
		
		
	}
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		//REAL(out)[i] = MAT_ELT(IX,i,1,n);
		REAL(out)[i] = IX[i][0];
	}
	for (j = 0; j < m; j++) {
		//REAL(out)[j+n] = MAT_ELT(IX,j,2,n);
		REAL(out)[j+n] = IX[j][1];
	}
	REAL(out)[m+n] = crt;
	return out;
	
}



// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//



// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//


SEXP barysort2(SEXP M, SEXP dims, SEXP pv , SEXP vs){
	
	int i,j;// rd, cd;//j2, i2, ki, kj
	
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	int nd = 2;
	
	//weight vectors
	float wr[n];
	float wc[m];
	
	for (i = 0; i < n; i++) {
		wr[i] = i;
	}
	for (i = 0; i < m; i++) {
		wc[i] = i;
	}
	
	//int MX[n][m];
	int nm = mymax(n,m);
	int IX[nm][2];
	int TX[nm][2];
	
	float rw[n];
	float cw[m];
	
	
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	int *mv = &INTEGER(M)[0];
   // int *IXp = &IX[0][0];
	int *TXp = &TX[0][0];
	
	//temporary row/col orders
	int tix[n];
	int tjx[m];
	
	
	
    //float* MX = calloc(n*m,sizeof(float));
    // one unnecessary copy here...
	
	float v[nm];
	//int no[nm];
	//int no2[nm];
	
	float* iv = &v[0];
	int* ptix = &tix[0];
	int* ptjx = &tjx[0];
	
	for (i = 0; i < n; i++) {
		rw[i] = 0;
	}
	
	for (j = 0; j < m; j++) {
		//MAT_ELT(IX,j,1,n) = j;
		//MAT_ELT(TX,j,1,n) = j;
		IX[j][1] = j;
		TX[j][1] = j;
		tjx[j] = j;
		cw[j] = 0;
		for (i = 0; i < n; i++) {
			//MAT_ELT(IX,i,2,n) = i;
			//MAT_ELT(TX,i,2,n) = i;
			IX[i][0] = i;
			TX[i][0] = i;
			tix[i] = i;
			//MAT_ELT(MX,i,j,n) = REAL(M)[i+j*n];//INTEGER
			cw[j] += mv[i+j*n];   //MAT_ELT(MX,i,j,n);
			rw[i] += mv[i+j*n]; //MAT_ELT(MX,i,j,n);
		}
	}
	int opt = 0;
	float crt0, crt1;
	float crt = 0;
	int row = 0;
	
	while(opt == 0){
		opt = 1;
		for (i=0; i < n; i++) {
			v[i] = 0;
			for (j = 0; j < m; j++) {
				//v[i] += MAT_ELT(MX,   MAT_ELT(IX,i,1,n), MAT_ELT(IX,j,2,n),n)*j;
				v[i] += mv[tix[i]+tjx[j]*n]*wc[j];   //MAT_ELT(MX,   tix[i], tjx[j],n)*j;
			}
			v[i] = v[i]/rw[tix[i]];
		}
		
		
		
		
		
		// get new row order
		quickSort(iv,ptix,0,n-1);
		
		//update row weights
		//for (i=1; i<n; i++) {
			//wr[i] = wr[i-1]+v[i]-v[i-1];
			
		//}
		for (i=0; i<n; i++) {
			wr[i] = v[i];
		}
		
		for (j=0; j < m; j++) {
			v[j] = 0;
			for (i = 0; i < n; i++) {
				v[j] += mv[tix[i]+tjx[j]*n]*wr[i];//MAT_ELT(MX,   tix[i], tjx[j],n)*i;
			}
			v[j] = v[j]/cw[tjx[j]];
		}
		
		// get new column order given the new row order
		quickSort(iv,ptjx,0,m-1);
		
		//update column weights
		//for (j=1; j<m; j++) {
		//	wc[j] = wr[j-1]+v[j]-v[j-1];
		//}
		for (j=0; j<m; j++) {
			wc[j] = v[j];
		}
		
		for (i = 0; i < n; i++) {
			TX[i][0] = tix[i];
		}
		// compute BCC after row order changes
		crt0 = mvclasscrit2( mv, TXp, dimv, nd);
		
		for (j = 0; j < m; j++) {
			TX[j][1] = tjx[j];
		}
		// compute BCC after row AND column order changes
		crt1 = mvclasscrit2( mv, TXp, dimv, nd);
		
		
		
		if (crt0 > crt) {
			opt = 0;
			crt = crt0;
			for (i=0; i<n; i++) {
				IX[i][0] = tix[i];
			}
			row = 1;
		}
		if (crt1 > crt) {
			opt = 0;
			crt = crt1;
			for (j=0; j<m; j++) {
				IX[j][1] = tjx[j];
			}
			row = 0;
		}
		
		//if (crt1 > crt | crt0 > crt) {
		//better than before... keep it and go on.
		//	opt = 0;
		//	crt = crt1;
		//	for (j=0; j<m; j++) {
		//MAT_ELT(IX,j,1,n) = MAT_ELT(TX,j,1,n);
		//MAT_ELT(IX,j,2,n) = MAT_ELT(TX,j,2,n);
		//		IX[j][1] = tjx[j];
		//	}
		//	for (i=0; i<n; i++) {
		//MAT_ELT(IX,j,1,n) = MAT_ELT(TX,j,1,n);
		//MAT_ELT(IX,j,2,n) = MAT_ELT(TX,j,2,n);
		//IX[i][0] = tix[i];
		//	}
		//}else {
		// do nothing and break;
		//}
		if( (crt1 < crt) & (crt0 < crt) ){
			// finish
			// weights are not correct in this part!!!
			if(row == 1){
				for (i = 0; i < n; i++) {
					tix[i] = IX[i][0];
				}
				for (j=0; j < m; j++) {
					v[j] = 0;
					for (i = 0; i < n; i++) {
						v[j] += mv[tix[i]+tjx[j]*n]*wr[i];//MAT_ELT(MX,   tix[i], tjx[j],n)*i;
					}
					v[j] = v[j]/cw[tjx[j]];
				}
				
				// get new column order given the new row order
				quickSort(iv,ptjx,0,m-1);
				for (j=0; j<m; j++) {
					IX[j][1] = tjx[j];
				}
			}else{
				// get best column order
				for (j = 0; j < m; j++) {
					tjx[j] = IX[j][1];
				}
				for (i=0; i < n; i++) {
					v[i] = 0;
					for (j = 0; j < m; j++) {
						//v[i] += MAT_ELT(MX,   MAT_ELT(IX,i,1,n), MAT_ELT(IX,j,2,n),n)*j;
						v[i] += mv[tix[i]+tjx[j]*n]*wc[j];   //MAT_ELT(MX,   tix[i], tjx[j],n)*j;
					}
					v[i] = v[i]/rw[tix[i]];
				}
				// reorder rows
				quickSort(iv,ptjx,0,m-1);
				for (i=0; i<n; i++) {
					IX[i][0] = tix[i];
				}
			}
		}
		
		
		
	}
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		//REAL(out)[i] = MAT_ELT(IX,i,1,n);
		REAL(out)[i] = IX[i][0];
	}
	for (j = 0; j < m; j++) {
		//REAL(out)[j+n] = MAT_ELT(IX,j,2,n);
		REAL(out)[j+n] = IX[j][1];
	}
	REAL(out)[m+n] = crt;
	return out;
	
}



// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//




// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//

	
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//



SEXP getclust(SEXP M, SEXP dims, SEXP tau0, SEXP method, SEXP singlesplit){
	
	int i,j,i2,j2, ki, kj, rd, cd;
	
	int v = INTEGER(method)[0];
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	
	
	//int MX[n][m];
    float* MX = calloc(n*m,sizeof(float));
    // one unnecessary copy here...
	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++) {
			MAT_ELT(MX,i,j,n) = REAL(M)[i+j*n];//INTEGER
		}
	}
	
	float currtau;
	float currbest;
	//long  a, b, c, d;
	float  a, b, c, d;
	float nn, zz;
	float p0, pc;
	int k = 0;
	int ncl = 1;
	int colcuts[m];
	int rowcuts[n];
	float tauval[n];
	int cutno[n];
	
	tauval[0] = 1;
	tauval[1] = 1;
	cutno[0] = 0;
	cutno[1] = 0;
	
	colcuts[0] = 0;
	rowcuts[0] = 0;
	colcuts[1] = m;
	rowcuts[1] = n;
	
	int cord[m];
	int rord[n];
	for (i=0; i<n; i++) {
		rord[i] = i;
	}
	for (j=0; j<m; j++) {
		cord[j]=j;
	}
	
	// diag prods for ties handling
	
	float cpr = 0;
		
	
	//int NE[n][m];
	//int NW[n][m];
	//int SE[n][m];
	//int SW[n][m];
    float* NE = calloc(n*m,sizeof(float));
    float* NW = calloc(n*m,sizeof(float));
    float* SE = calloc(n*m,sizeof(float));
    float* SW = calloc(n*m,sizeof(float));
	
	for (j = 0; j < m; j++) {
		MAT_ELT(NE,0,j,n) = MAT_ELT(MX,0,j,n);
		MAT_ELT(SE,n-1,j,n) = MAT_ELT(MX,n-1,j,n);
		MAT_ELT(NW,0,j,n) = MAT_ELT(MX,0,j,n);
		MAT_ELT(SW,n-1,j,n) = MAT_ELT(MX,n-1,j,n);
	}

	for (i = 1; i < n; i++) {
		for (j = 0; j < m; j++) {
			MAT_ELT(NW,i,j,n) = MAT_ELT(NW,i-1,j,n) + MAT_ELT(MX,i,j,n);
			MAT_ELT(NE,i,j,n) = MAT_ELT(NW,i-1,j,n) + MAT_ELT(MX,i,j,n);
			MAT_ELT(SE,n-i-1,j,n) = MAT_ELT(SE,n-i,j,n) + MAT_ELT(MX,n-i-1,j,n);
			MAT_ELT(SW,n-i-1,j,n) = MAT_ELT(SW,n-i,j,n) + MAT_ELT(MX,n-i-1,j,n);
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 1; j < m; j++) {
			MAT_ELT(NW,i,j,n) = MAT_ELT(NW,i,j-1,n) + MAT_ELT(NW,i,j,n);
            MAT_ELT(NE,i,m-j-1,n) = MAT_ELT(NE,i,m-j,n) + MAT_ELT(NE,i,m-j-1,n);
            MAT_ELT(SE,i,m-j-1,n) = MAT_ELT(SE,i,m-j,n) + MAT_ELT(SE,i,m-j-1,n);
            MAT_ELT(SW,i,j,n) = MAT_ELT(SW,i,j-1,n) + MAT_ELT(SW,i,j,n);
		}
	}

	
	
	///////////
	// a | b //
	//-------//
	// c | d //
	///////////
	float rs1,rs2,cs1,cs2;
	while ( (k < ncl) ) {
		//currbest = -0.001;
		cpr = 0;
		currbest = -1.0001;
		rd = rowcuts[k+1]-rowcuts[k];
		cd = colcuts[k+1]-colcuts[k];
		
		if( (rd > 1)  && ( cd > 1)   ){
			
			for (i2 = rowcuts[k]+1; i2 < rowcuts[k+1]; i2++) { // i2 and j2 are first in the second part
				for (j2 = colcuts[k]+1; j2 < colcuts[k+1]; j2++) {
					
					a = MAT_ELT(NW,i2-1,j2-1,n);
					b = MAT_ELT(NE,i2-1,j2,n);
					c = MAT_ELT(SW,i2,j2-1,n);
					d = MAT_ELT(SE,i2,j2,n);
					if( k > 0 ){//if (colcuts[k] > 0) { //left+top
						a = MAT_ELT(NW,i2-1,j2-1,n) - MAT_ELT(NW,i2-1, colcuts[k]-1 ,n) - MAT_ELT(NW, rowcuts[k]-1 ,j2-1,n) + MAT_ELT(NW, rowcuts[k]-1 , colcuts[k]-1 ,n);
						c = c - MAT_ELT(SW,i2, colcuts[k]-1 ,n);
						b = b - MAT_ELT(NE, rowcuts[k]-1 , j2 ,n);
					}
					if( k+1 < ncl ){//if (colcuts[k+1] < m-1) { //right+bottom
						d = MAT_ELT(SE,i2,j2,n) - MAT_ELT(SE,i2, colcuts[k+1] ,n) - MAT_ELT(SE, rowcuts[k+1] ,j2,n) + MAT_ELT(SE, rowcuts[k+1] , colcuts[k+1] ,n);
						c = c - MAT_ELT(SW, rowcuts[k+1] ,j2-1,n); 
						b = b - MAT_ELT(NE, i2-1 , colcuts[k+1] ,n);
					}
					if((k > 0) && (k+1 < ncl)){//if (colcuts[k+1] < m-1 && colcuts[k] > 0) { //all
						c = c + MAT_ELT(SW, rowcuts[k+1] , colcuts[k]-1 ,n);
						b = b + MAT_ELT(NE, rowcuts[k]-1 , colcuts[k+1] ,n);
					}
					//Rprintf(" a = %d, b = %d, c = %d, d = %d\n",a,b,c,d);
				if(v == 1){	
					//kendalls
					zz = a*d - b*c;
					nn = pow((a+b)*(c+d),0.5)*pow((a+c)*(b+d),0.5);
				}
				if(v == 2){	
					//kappa
					pc = (float)((a+b)*(a+c)+(c+d)*(b+d))/(a+b+c+d)/(a+b+c+d);
					p0 = (float)(a+d)/(a+b+c+d);
				  zz = p0-pc;
				  nn = 1-pc;
				}
					if(v == 3){	
						//WBCI
						rs1 = a+b;
						rs2 = c+d;
						cs1 = a+c;
						cs2 = b+d;
						pc = (float)(c*(a+d+2*b) + a*b + d*b);
						p0 = (float)( rs2*cs1*( rs1*cs1 +rs2*cs2 + 2*rs1*cs2  ) + rs1*cs1*rs1*cs2 + rs2*cs2*rs1*cs2  )/(rs1+rs2)/(rs1+rs2);
						zz = p0 - pc;
						nn = p0;
					}	
				
					if(v == 4){	
						//BCI
						nn = (a+b)*(c+d)*(a+c)*(b+d);
						nn = nn / ( (a+b+c+d)*(a+b+c+d) );
						zz = nn - b*c;
					}
					if(v == 5){	
						rs1 = a+b;
						rs2 = c+d;
						cs1 = a+c;
						cs2 = b+d;
						nn = a+b+c+d;
						p0 = (nn*c - rs2*cs1)/pow(rs2*cs1,0.5);
						pc = (nn*b - rs1*cs2)/pow(rs1*cs2,0.5);
						zz =  mymin2( p0, pc);
						zz = log(2) - log(1 + exp(zz/pow(nn,0.5)/4));
						nn = log(2);
					}
                    if(v == 6){
                        //Rmin
						nn = a+b+c+d;
						rs1 = (a+b)/nn;
						rs2 = (c+d)/nn;
						cs1 = (a+c)/nn;
						cs2 = (b+d)/nn;
                        
                        p0 = (a/nn - rs1*cs1)/pow(rs1*cs1,0.5);
						
                        p0 = mymin2( p0 , (b/nn - rs1*cs2)/pow(rs1*cs2,0.5) );
						p0 = mymin2( p0 , (c/nn - rs2*cs1)/pow(rs2*cs1,0.5) );
						p0 = mymin2( p0 , (d/nn - rs2*cs2)/pow(rs2*cs2,0.5) );
						
						
						zz =  -p0;
						nn = 1;
					}
				  // Rprintf(" zz = %f, nn = %f\n",zz,nn);
					
					currtau = zz/nn;
				
					if( (currtau == currbest) & ( a*d > cpr) ){
						// prefer balanced cuts
						//Rprintf("prefering %d, %d over %d, %d with cellprod = %f",i2,j2,ki,kj,cpr);
						cpr = a*d;
						currbest = currtau;
						ki = i2;
						kj = j2;
					}
					
					
					if(currtau > currbest){
						currbest = currtau;
						ki = i2;
						kj = j2;
						
					}
				}
			}
			//Rprintf("best crit = %f at x = %d and y = %d\n",currbest,kj,ki);
			
			if(currbest >= REAL(tau0)[0]){
				ncl++;
			
				for (i = ncl; i > k+1; i--) {
					colcuts[i] = colcuts[i-1];
					tauval[i] = tauval[i-1];
					cutno[i] = cutno[i-1];
				}
				colcuts[k+1] = kj;
				tauval[k+1] = currbest;
				cutno[k+1] = mymax(cutno[k+2],cutno[k])+1;
				for (i = ncl; i > k+1; i--) {
					rowcuts[i] = rowcuts[i-1];
				}
				rowcuts[k+1] = ki;
				if(ncl == INTEGER(singlesplit)[0]+1){
					k = ncl;
				}
			}else {
				k++;
			}
			
		}else {
			k++;
		}

		
				
		
	}
	
	SEXP out = allocVector(REALSXP,4*k);
	for( i = 0; i < k;i++ ){
		REAL(out)[i] = rowcuts[i+1];
		REAL(out)[i+k] = colcuts[i+1];
		REAL(out)[i+2*k] = tauval[i+1];
		REAL(out)[i+3*k] = cutno[i+1];
	}
    free(NE);
    free(NW);
    free(SE);
    free(SW);
    free(MX);
	return out;
	
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //





// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //



float mvhammcrit(int *m0, int *IM, int dims[], int nd){
	
	
	int i, j1,j2, s1, s2, k, d, rind ,lind;//s, j

	int (*IX)[nd] = (int (*)[nd])IM;
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;

	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dims[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dims[nd-i];//CHECK
	}
	// m1 array for the cumsums
	int ml = cp[nd-1]*dims[nd-1];
	int m1[ ml ];
	int m2[ ml ];
	//int m2[ ml ];
	for (i = 0; i < ml; i++) {
		m1[i] = m0[i];
		m2[i] = 0;
	}
	//Rprintf("\n\n ml = %d\n\n",ml);

	//go through m1 and m0 with stepsize s acc. to dimension k
	// dims is in reverse order
	
	for (d = 0; d < nd; d++) {
		
			s1 = cp[d];
			s2 = cpr[d];
			
			// first k-vector is zero
			for (j1 = 0; j1 < s2; j1++) {
				for (j2 = 0; j2 < s1; j2++) {
					rind = IX[0][d]*s1 + j1*s1*dims[d] + j2;
					m1[ rind ] = 0;
				}
			}
			
			// cumsum over k = d 1st
			for (i = 1; i < dims[d]; i++) {
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						rind = IX[i][d]*s1 + j1*s1*dims[d] + j2;
						lind = IX[i-1][d]*s1 + j1*s1*dims[d] + j2;
						m1[ rind ] = m1[ lind ] + m0[ lind ];
					}
				}
			}
		
			// cumsum over k = d 2nd
			for (i = 1; i < dims[d]; i++) {
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						rind = IX[i][d]*s1 + j1*s1*dims[d] + j2;
						lind = IX[i-1][d]*s1 + j1*s1*dims[d] + j2;
						m1[ rind ] += m1[ lind ];
					}
				}
			}
		
		// cumsum over other dimensions once
		// add to m2 and reset m1
		for (k = 0; k < nd; k++) {
			s1 = cp[k];
			s2 = cpr[k];
			if (d != k) {
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						rind = IX[0][k]*s1 + j1*s1*dims[k] + j2;
						m2[ rind ] += m1[ rind ];
						for (i = 1; i < dims[k]; i++) {
							rind = IX[i][k]*s1 + j1*s1*dims[k] + j2;
							lind = IX[i-1][k]*s1 + j1*s1*dims[k] + j2;
							m1[ rind ] += m1[ lind ];
							m2[ rind ] += m1[ rind ];
						}
					}
				}
			} 
		}
		for (i = 0; i < ml; i++) {
			m1[i] = m0[i];
		}
	}
	float cv = 0.0;
	for (i = 0; i < ml; i++) {
		cv += m2[i] * m0[i];
	}
	return cv;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP classcrit(SEXP M, SEXP dims, SEXP vs){
	int i,k;
	int nd = LENGTH(dims);
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	
	int mt = INTEGER(vs)[0];
	int* mv = &INTEGER(M)[0];
	
	int ml = dimv[0];
	for(i = 1; i < nd; i++){
		if( ml < dimv[i] )
			ml = dimv[i];
	}
	
	int IM[ml][nd];
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IM[i][k] = i;
		}
	}
	int* IMp = &IM[0][0];
	float crit;
	if(mt == 0){
		crit = mvclasscrit2(mv, IMp, dimv, nd);
	}else{
		crit = mvhammcrit(mv, IMp, dimv, nd);
	}
	SEXP out = allocVector(REALSXP, 1);
	REAL(out)[0] = crit;
	return out;
	
}


// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //



// this function will perform permutations for the categories until a local optimum is reached
SEXP mvclass(SEXP M, SEXP dims, SEXP pv , SEXP vs){
	
	
	int s,i,j,k,ii,jj,kk, tmp;//r, t
	int nd = LENGTH(dims);
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	int mt = INTEGER(vs)[0];
	
	int ml = dimv[0];
	int biglen = ml;
	for(i = 1; i < nd; i++){
		biglen = biglen + dimv[i];
		if( ml < dimv[i] )
			ml = dimv[i];
	}
	int IM[ml][nd];
	int TIM[ml][nd];
	
	//Rprintf("init done");
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IM[i][k] = i;
			TIM[i][k] = i;
			//Rprintf("%d\t",TIM[i][k]);
		}
		//Rprintf("\n");
	}
	
	
	
	int* mv = &INTEGER(M)[0];
	int* TIMp = &TIM[0][0];
	
	
	int bestmove[3];
	float bestcrit = 0;//FLT_MAX;
	float tempcrit;
	int better;
	int opt = 0;
	//Rprintf("bestcrit = %f\n",bestcrit);
	
	while (opt == 0) {
		
		better = 0;
		for (k = 0; k < nd; k++) {
			
			
			if (INTEGER(pv)[k] == 1) {
				
				for (i = 0; i < dimv[k]; i++) { // move i to j in dim k

					for (j = 0; j < dimv[k]; j++) {
						//Rprintf("k = %d\t i = %d \tj = %d\t",k,i,j);
						if (i < j) {
							for (s = i; s < j; s++) {
								TIM[s][k] = IM[s+1][k];
							}
							TIM[j][k] = IM[i][k];
						}
						if (i > j) {
							for (s = i; s > j; s--) {
								TIM[s][k] = IM[s-1][k];
							}
							TIM[j][k] = IM[i][k];
						}		
					
						if( mt == 0 ){
							tempcrit =  mvclasscrit2(mv, TIMp, dimv, nd);
						}
						if( mt == 1 ){
							tempcrit =  mvhammcrit(mv, TIMp, dimv, nd);
						}
												
						//Rprintf("tempcrit = %f\n", tempcrit);
						if (tempcrit > bestcrit) {
							bestcrit = tempcrit;
							bestmove[0] = k;
							bestmove[1] = i;
							bestmove[2] = j;
							better = 1;
						//	Rprintf("better = %f -----at %d  %d  %d------------ >>><<<\n", bestcrit,k,i,j);
						}
						//reset TIM
						if (i < j) {
							for (s = j; s >= i; s--) {
								TIM[s][k] = IM[s][k];
							}
						}
						if (i > j) {
							for (s = i; s >= j; s--) {
								TIM[s][k] = IM[s][k];
							}
						}

						
						
						
					} // j
				} // i
			} // if
		} // k
		
		if(better > 0){
			// changes to IM
			ii = bestmove[1];
			jj = bestmove[2];
			kk = bestmove[0];
			//Rprintf("permutation dim = %d, %d to %d\n ------->>><<<----->>><<<----#+###+#",kk,ii,jj);
			tmp = IM[ii][kk];
			if (ii < jj) {
				for (s = ii; s < jj; s++) {
					IM[s][kk] = IM[s+1][kk];
					TIM[s][kk] = IM[s][kk];
				}
				IM[jj][kk] = tmp;
				TIM[jj][kk] = IM[jj][kk];
			}
			if (ii > jj) {
				for (s = ii; s > jj; s--) {
					IM[s][kk] = IM[s-1][kk];
					TIM[s][kk] = IM[s][kk];
				}
				IM[jj][kk] = tmp;
				TIM[jj][kk] = IM[jj][kk];
			}
			
			//Rprintf("\n The new orders are:\n");
			//for (r = 0; r < nd; r++) {
			//	for (t = 0; t < dimv[r]; t++) {
			//		Rprintf("%d, ",IM[t][r]);
			//	}
			//	Rprintf("\n");
			//}
			
			
		}else{
			opt = 1;	
		}	
	}
	
	SEXP out = allocVector(REALSXP, biglen + 1);
	
	k = 0;
	for(i = 0; i < nd; i++){
		for (j = 0; j < dimv[i]; j++) {
			REAL(out)[k] = IM[j][i];
			k++;
		}
	}
	REAL(out)[biglen] = bestcrit;
	
	return out;
}










// this function will perform permutations for the categories until a local optimum is reached
SEXP preclass(SEXP M, SEXP dims, SEXP pv , SEXP vs, SEXP sym){
	
	
	//Rprintf("SYM = %d\n",INTEGER(sym)[0]);
	
	int s,i,j,k, d;//r, tmp, t, kk, ii, jj
	int nd = LENGTH(dims);
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	
	
	int ml = dimv[0];
	int biglen = ml;
	for(i = 1; i < nd; i++){
		biglen = biglen + dimv[i];
		if( ml < dimv[i] )
			ml = dimv[i];
	}
	int IM[ml][nd];
	int CIM[ml][nd]; // current IM
	int TIM[ml][nd]; // a dummy order matrix
	
	int* CIMp = &CIM[0][0];
	int* TIMp = &TIM[0][0];
	int* IX = &IM[0][0];
	
	float cvv[ml];
	float *cvp = &cvv[0];
	
	//Rprintf("init done");
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IM[i][k] = i;
			CIM[i][k] = i;
		}
	}
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < ml; i++) {
			TIM[i][k] = i;
		}
	}
	
	int* mv = &INTEGER(M)[0];
		
	float bestcrit = 0;//FLT_MAX;
	float newcrit, oldcrit;
	
	int opt = 0;
	//Rprintf("bestcrit = %f\n",bestcrit);
	
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dimv[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dimv[nd-i];//CHECK
	}
	// m1 array for the cumsums
	int mln = cp[nd-1]*dimv[nd-1];
	//int m1[ mln ];
	//int m1i[ mln ];
	
	//int m0[ mln ];
	int* m1 = calloc(mln,sizeof(int));
	int* m1i = calloc(mln,sizeof(int));
	int* m0 = calloc(mln,sizeof(int));
	
	for (i = 0; i < mln; i++) {
		m0[i] = INTEGER(M)[i];
	}
	
	for (i = 0; i < mln; i++) {
		m1[i] = m0[i];
		m1i[i] = 0;
	}
	
	//int* m1ip = &m1i[0];
	
	int ndt = nd;
	if(INTEGER(sym)[0] == 1){
		ndt = 1;
	}
	
	int ds, s1, s2, j1, j2, rind, lind;//s1t, s2t
	int tmpdimv[nd];
	int iweight;
	int ord[ml];
	int *ordp = &ord[0];
	
	while (opt == 0) {
		//////Rprintf("hoi\n\n");
		opt = 1;
		for (d = 0; d < ndt; d++) {
			
			
			if (INTEGER(pv)[d] == 1) {
					
					s1 = cp[d];
					s2 = cpr[d];
					
					// first k-vector is zero
					for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
							rind = IM[0][d]*s1 + j1*s1*dimv[d] + j2;
							m1[ rind ] = m0[ rind ];
						}
					}
								
				
					// sum over k = d => the full sums are in i = dims[d]-1
					for (i = 1; i < dimv[d]; i++) {
						for (j1 = 0; j1 < s2; j1++) {
							for (j2 = 0; j2 < s1; j2++) {
								rind = IM[i][d]*s1 + j1*s1*dimv[d] + j2;
								lind = IM[i-1][d]*s1 + j1*s1*dimv[d] + j2;
								m1[ rind ] = m1[ lind ] + m0[ rind ];
							}
						}
					}
				
				
				//writing the sums to m1i
				rind = 0;
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						//rind = 0*s1 + j1*s1*dimv[d] + j2;
						lind = IM[dimv[d]-1][d]*s1 + j1*s1*dimv[d] + j2;
						m1i[ rind*2 ] = m1[ lind ];
						//Rprintf("val %d from %d to %d\n", m1[ lind ],lind,rind);
						rind++;
					}
				}
				
				//temporary dimv
				ds = 1;
				tmpdimv[0] = 2;
				
				for (s = 0; s < 2; s++) {
					TIM[s][0] = s;
				}
				for (k = 0; k < nd; k++) {
						if (k != d) {
							for (s = 0; s < dimv[k]; s++) {
								TIM[s][ds] = IM[s][k];
							}
							tmpdimv[ds] = dimv[k];
							ds++;
						}
				}
							
				s1 = cp[d];
				s2 = cpr[d];
				
				for (i = 0; i < dimv[d]; i++) {
					rind = 0;
					// writing category i to m1i
					iweight = 0;
					for (j2 = 0; j2 < s1; j2++) {
						for (j1 = 0; j1 < s2; j1++) {
							//IS THIS CORRECT??? IX[i][d] or simply i ????
							lind = IM[i][d]*s1 + j1*s1*dimv[d] + j2;
                            //lind = i*s1 + j1*s1*dimv[d] + j2;
							m1i[ rind*2+1 ] = m0[ lind ];
							iweight += m0[ lind ];
                            rind++;
						}
					}
			
										
					cvv[i] = mvclasscrit2(m1i, TIMp, tmpdimv, nd) / iweight;//here we had the old mvclasscrit and m1ip
					//Rprintf("cvv[%d] = %f\n",i,cvv[i]);
				}	
				
				
			}
			for (s = 0; s < dimv[d]; s++) {
				ord[s] = CIM[s][d];
			}
			
			quickSort(cvp,ordp,0,dimv[d]-1);	
			
            //CIM is the temporary order resulting from the quicksort
			for (s = 0; s < dimv[d]; s++) {
				//TEST: odd dimensions (1,3,...) inverse??
				//if( (d % 2) == 1){
				//	CIM[s][d] = ord[dimv[d]-s-1];
				//}else {
					CIM[s][d] = ord[s];
				//}

			}
			
						
			oldcrit = mvclasscrit2(mv, IX, dimv, nd);//here we had the old mvclasscrit
			//Rprintf("\n\ndim = %d, old crit = %f\n\n", d,oldcrit);
			//Rprintf("c(");
			//for (i = 0; i < dimv[d]-1; i++) {
			//	Rprintf("%d ,",CIM[i][d]);
			//}
			//Rprintf("%d )\n\n",CIM[dimv[d]-1][d]);
			
			if(INTEGER(sym)[0] == 1){
				for (j = 0; j < nd; j++) {
					for (i = 0; i < dimv[j]; i++) {
						CIM[i][j] = ord[i];
					}
				}
			}
			
			newcrit = mvclasscrit2(mv, CIMp, dimv, nd);
			
			if (newcrit > bestcrit) {
				bestcrit = newcrit;
				//////Rprintf("new best crit = %f", bestcrit);
				
				if(INTEGER(sym)[0] == 1){
					for (j = 0; j < nd; j++) {
						for (i = 0; i < dimv[j]; i++) {
							IM[i][j] = CIM[i][j];
						}
					}
				}else{
					for (s = 0; s < dimv[d]; s++) {
						IM[s][d] = CIM[s][d];
					}
				}
				opt = 0;
			}else {
				//////Rprintf("newcrit = %f stop", newcrit);
				// ????
				if(INTEGER(sym)[0] == 1){
					for (j = 0; j < nd; j++) {
						for (i = 0; i < dimv[j]; i++) {
							CIM[i][j] = IM[i][j];
						}
					}
				}else{
					for (s = 0; s < dimv[d]; s++) {
						CIM[s][d] = IM[s][d];
					}
				}
				
			}
		}
	}
				
	SEXP out = allocVector(REALSXP, biglen + 1);
	
	k = 0;
	for(i = 0; i < nd; i++){
		for (j = 0; j < dimv[i]; j++) {
			REAL(out)[k] = IM[j][i];
			k++;
		}
	}
	REAL(out)[biglen] = bestcrit;
	free(m1);
    free(m0);
    free(m1i);
	return out;
}



////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////







SEXP cumsum(SEXP M, SEXP dimv){
	
		int i, j1,j2,s, s1, s2, k,lind, d, ind1, ind2;//ind3, tmp, tmp2, j, rind, ix
	
	int nd = LENGTH(dimv);
	
	int dims[nd];
	for (i=0; i< nd; i++) {
		dims[i] = INTEGER(dimv)[i];
	}

	int *m0 = &INTEGER(M)[0];
	
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dims[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dims[nd-i];//CHECK
	}
	
	// m1 array for the cumsums
	int ml = cp[nd-1]*dims[nd-1];
	int m1[ ml ];
	int m2[ ml ];
	//int m2[ ml ];
	for (i = 0; i < ml; i++) {
		m1[i] = m0[i];
		//m2[i] = 0;
	}
	

		
		// cumsum over all dimensions with a shift of 1:
    // 2new = 1old + 0new, 3new = 2old + 1new, ... => save (e.g.) 2old for i = 3 step in tmp(2)
	
// compute the cumulative sums		
 for (k = 0; k < nd; k++) {
            s1 = cp[k];
            s2 = cpr[k];
                for (j1 = 0; j1 < s2; j1++) {
                    for (j2 = 0; j2 < s1; j2++) {
                        for (i = 1; i < dims[k]; i++) {
                            ind2 = i*s1 + j1*s1*dims[k] + j2;
                            ind1 = ind2 - s1;
                            m1[ ind2 ] += m1[ ind1 ];
                      	}
                    }
                }
            
		} 
// write the cumulative sums to index+1 in m2
	// compute the increment in m0 and m1
	int cpx[nd];
	cpx[0] = dims[0];
	
	for (i = 1; i < nd; i++) {
		cpx[i] = cpx[i-1]*dims[i];
   }
	
    
	int inc = 1;
	for (s = 0; s < nd-1; s++) {
		inc += cpx[s];
	}
   // Rprintf("%d\n",inc);
	for (i = 0; i < ml-inc; i++) {
        m2[i + inc] = m1[i];
    }

    // set the border (any dim = 0)
	for (d = 0; d < nd; d++) {
		s1 = cp[d];
		s2 = cpr[d];
		for (j1 = 0; j1 < s2; j1++) {
			for (j2 = 0; j2 < s1; j2++) {
				lind = 0*s1 + j1*s1*dims[d] + j2;
               m2[ lind ] = 0;
			}
		}
	}

	
	SEXP out = allocVector(REALSXP, ml);
	for (i = 0; i < ml; i++) {
		REAL(out)[i] = m2[i];
	}
	
	return out;
}




// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////







// this function will perform permutations for the categories until a local optimum is reached
SEXP quickhamm(SEXP M, SEXP dims, SEXP pv , SEXP vs, SEXP minc, SEXP treevec, SEXP treelengths){
	
	// the input is as follows:
	
	// > M is the array which is to be optimized
	// > dims is a vector with the numbers of categories for each variable
	// > pv is a 0/1 vector which indicates whether or not a variable shall be reordered
	// > vs is a (currently unused) version number
	// > tree contains a (binary) tree object with 4 entries for each node:
	//  the first and the last index of the left and the right branch
	//  usually the right index of the left branch is 1 less than the left index of the right branch
	// > CM is a nd x nd matrix with optimal values for the pairwise criteria. It is used to decide which variable shall be optimized in each major step.
	// > minc (eps) is the minimal improvement which is required for a movement as a ratio of the total criterion value at the beginning of a major step.
	//  i.e. if eps = 0.01 only movements with an improvement of at least 1 percent are possible.
	
	int s,t,i,j,k,tmp, tmp2, d, g0, g1, g2, i1, i2, g;//r
	int nd = LENGTH(dims);
    
    tmp = INTEGER(treelengths)[0];
    
	int tl[nd];
    for (i = 0; i < nd; i++) {
        if (LENGTH(treelengths) > 1 ) {
            tl[i] = INTEGER(treelengths)[i];
        }else{
            tl[i] = 0;
        }
    }
    
    
    if( LENGTH(treelengths) > 1 ){	
        // set tree parameters and objects
        for (i = 1; i < nd; i++) {
            tmp = mymax(tmp,INTEGER(treelengths)[i]);
        }
        tmp = 3*tmp;
    }
    //int TM[tmp][nd];
    int *TM = calloc(tmp*nd,sizeof(int));
    
    //int *tl = &INTEGER(treelengths)[0];
	// TM is equivalent to trees, but easier to read in the code later
	// it contains the indices of the node children: min.left, min.right = max.left +1, max.right
	// TM will not change during the algorithm, the changes are made in the index vector it refers to
    
    if( LENGTH(treelengths) > 1 ){	
        k = 0;
        for (i = 0; i < tmp; i++) {
            for (j = 0; j < nd; j++) {
                MAT_ELT(TM,i,j,tmp) = 0;//TM[i][j] = 0;
            }
        }
        for (d = 0; d < nd; d++) {
            for (i = 0; i < tl[d]; i++) {
                MAT_ELT(TM,3*i,d,tmp) = INTEGER(treevec)[k + 3*i]-1;
                MAT_ELT(TM,3*i+1,d,tmp) = INTEGER(treevec)[k + 3*i+1]-1;
                MAT_ELT(TM,3*i+2,d,tmp) = INTEGER(treevec)[k + 3*i+2]-1;
            }
            k += 3*tl[d];
        }
    }
    
    // Rprintf("TM = \n");
	//for (i = 0; i < tmp; i++) {
	//	for (j = 0; j < nd; j++) {
	//		Rprintf("%d\t",TM[i][j]);
	//	}
    //   Rprintf("\n");
	//}
    
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	//Rprintf("no. of trees = %d", nt);  
	// ml is the highest number of categories among all variables
	// biglen is the total number of categories ober all variables
	int ml = dimv[0];
	int biglen = ml;
	for(i = 1; i < nd; i++){
		biglen = biglen + dimv[i];
		if( ml < dimv[i] )
			ml = dimv[i];
	}
    
	// IX is a matrix with the current order of indices. Each column is for one variable
	// PX contains the current position of the original indices (only for trees)
	int IX[ml][nd];
	int PX[ml][nd];
    int cind[ml];
	// initial indices
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IX[i][k] = i;
			PX[i][k] = i;
		}
	}
    
    int *mv = &INTEGER(M)[0];
    int *IXp = &IX[0][0];
	// cp and cpr are products of the numbers of categories for the variables before and after dimension d
	// this is needed to grab the right entries from the m0 and m1 vectors
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dimv[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dimv[nd-i];
	}
	
	// the m1 array will contain the temporary (for one major step) cumulative sums over all dimensions except d
	// m0 contains the original values and is used to feed m1 in each step
	// m2 is the same as m1 except that all entries are relocated at one position further (i -> i+1)
	// all three have length mln
	int mln = cp[nd-1]*dimv[nd-1];
	int m1[ mln ];
	int m0[ mln ];
	int m2[ mln ];
	
	for (i = 0; i < mln; i++) {
		m0[i] = INTEGER(M)[i];
	}
	for (i = 0; i < mln; i++) {
		m1[i] = m0[i];
		m2[i] = 0;
	}
    
	// TODO: compute the overall criterion alongside the computations
	// (init step and inc)
	
	float bestcrit = 0;//FLT_MAX;
	float newcrit = 0;
	//float startcrit = 0;
	newcrit = mvclasscrit2(mv, IXp, dimv, nd);
    //Rprintf("initial crit = %f",newcrit);
	
	int opt = 0;
	
	int  s1, s2,  j1, j2, rind, lind;//ds, s1t, s2t,
    
	// DM contains the pairwise delta values
	// (i.e. the difference between the criterion when cat. 1 is before or after cat. 2)
	// DMC contains the cumulative sums over rows of DM starting from the diagonal outwards
	int DM[ml][ml];
	long DMC[ml][ml];	 
    long DMCw[ml][ml];	
    long DMCl[ml][ml];	
    
    // best is the best movement delta value from DMC
    // ii and jj are the corresp. indices
    // eps is the minimal value for a movement
    // eps should be either 0 or decrease with every major step
    
    float best;
    int ii, jj,  ll, lr;//gg
    float eps = 0;
	float imp;
    
	/// init done
	//Rprintf("begin optim\n");
	// d is the current dimension considered for optimization
	d = 0;
	int init = 1;
    
    // ___   MAIN STEP   ___ //
    //  |				  |  //
    //  V                 V  //   
    
	while(opt == 0){
		// optimality is assumed unless an improvement can be made
		opt = 1;
		
        // current index for this dimension to work with
        for (i = 0; i < ml; i++) {
            cind[i] = i;
        }
        for (i = 0; i < ml; i++) {
            for (j = 0; j < ml; j++) {
                DM[i][j] = 0;
                DMC[i][j] = 0;
                DMCw[i][j] = 0;
                DMCl[i][j] = 0;
            }
        }
		// cumsum over all dimensions except d
        for (k = 0; k < nd; k++) {
			s1 = cp[k];
			s2 = cpr[k];
			if (d != k) {
				for (j1 = 0; j1 < s2; j1++) {
                    for (j2 = 0; j2 < s1; j2++) {
                        for (i = 1; i < dimv[k]; i++) {
                            //rind = IX[i][k]*s1 + j1*s1*dimv[k] + j2;
                            //lind = IX[i-1][k]*s1 + j1*s1*dimv[k] + j2;
                            rind = i*s1 + j1*s1*dimv[k] + j2;
                            lind = rind-s1;
							m1[ rind ] += m1[ lind ];
						}
						
					}
				}
			}
		} 
        
        
        //     Rprintf("test1");
        // shift the cumsums to index+1 (except dim d)
        // >> compute index increment
        int cpx[nd];
        cpx[0] = dimv[0];
        
        for (i = 1; i < nd; i++) {
            cpx[i] = cpx[i-1]*dimv[i];
        }
        int inc;
        if(d == 0){
            inc = 0;
        }else{
            inc = 1;
        }
        for (s = 0; s < nd-1; s++) {
            if( s != d-1 ){
                inc += cpx[s];
            }
        }
        //   Rprintf("inc = %d\n",inc);
        for (i = 0; i < mln-inc; i++) {
            m2[i + inc] = m1[i];
        }
        //  Rprintf("test2");
        // >> set the border (any dim expet d) to 0
        for (k = 0; k < nd; k++) {
            if(k != d){
                s1 = cp[k];
                s2 = cpr[k];
                for (j1 = 0; j1 < s2; j1++) {
                    for (j2 = 0; j2 < s1; j2++) {
                        //lind = IX[0][k]*s1 + j1*s1*dimv[k] + j2;
                        lind =  j1*s1*dimv[k] + j2;
                        m2[ lind ] = 0;
                    }
                }
            }
        }
        //  Rprintf("cumsums done");
        // for (i = 0; i < mln; i++) {
        // Rprintf(" %d ",m2[i]);
        //} 
		// compute pairwise criteria from m2 and m0
		// the mv criterion is to maximize the number of 'diagonal' pairs
		// i.e. the cumsums 
        // >> first write 'i1 is before i2 to DM
        // >> (we lose that much if we put i1 after i2)
        
        s1 = cp[d];
        s2 = cpr[d];
		for (i1 = 0; i1 < dimv[d]; i1++) {
			for (i2 = 0; i2 < dimv[d]; i2++) {
				if(i1 != i2){
					for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
							//lind = IX[i1][d]*s1 + j1*s1*dimv[d] + j2;
							//rind = IX[i2][d]*s1 + j1*s1*dimv[d] + j2;
                            lind = i1*s1 + j1*s1*dimv[d] + j2;
							rind = i2*s1 + j1*s1*dimv[d] + j2;
							DM[i1][i2] += m0[ rind ] * m2[ lind ]; 
						}
					}
				}
			}
		}
    
        // >> now write 'put i1 after i2 difference' to DM
		
        
        //for (i1 = 0; i1 < dimv[d]-1; i1++) {
		//	for (i2 = i1; i2 < dimv[d]; i2++) {
		//		if( (init == 1) & (i2 > i1) ){
		//			startcrit += DM[i1][i2];
		//		}
         //       DM[i1][i2] = DM[i2][i1] - DM[i1][i2];
         //       DM[i2][i1] = -DM[i1][i2];
        //    }
		//}
		//if( (init == 1) ){
		//	Rprintf("startcrit = %f\n", startcrit);
			init = 0;
		//}
		
    
        
		// reset m1 to m0,(m2 to 0?)
		for (i = 0; i < mln; i++) {
            m1[i] = m0[i];
            m2[i] = 0;
		}
        
        
		
        /// -----------------------------------------------------------------------		
        
		// ___   INTRA STEP 2  ___ //
		//  |				    |  //
		//  V                   V  //  
		
		// compute optimal order for veriable d from DM
		// note that DM is in the original order, DMC is in IX order
        
		// improvement must be better than eps to proceed
        imp = eps + 1;
		while( imp > eps ){
			
			if(tl[d] > 0){
                // ---------------- TREE STEP ---------------- //
				// there is a tree for dimension d	
				//Rprintf("treestep\n\n");
				// use IX = PX.old (IX is overwritten at the end of the tree step)
                // PX.length = dimv or tl?
             //   Rprintf("tree %d", d);
                
				for (i = 0; i < dimv[d]; i++) {
                    cind[i] = IX[i][d];
					IX[i][d] = PX[i][d];
				}
				//anybetter = 0;
                
				// go through all nodes
				// both DM and the tree work on the original indices
                // i.e. we always use the same indices, but on a changed DM matrix
                // save the differences in DMC[0][ 0..tl[d]-1 ]
                for (g = 0; g < tl[d]; g++) {
					DMC[ 0 ][ g ] = 0;
					g0 = MAT_ELT(TM,3*g,d,tmp);//TM[ g*3 + 0][d];
					g1 = MAT_ELT(TM,3*g+1,d,tmp);//TM[ g*3 + 1][d];
					g2 = MAT_ELT(TM,3*g+2,d,tmp);//TM[ g*3 + 2][d];
					//Rprintf("check node (%d, %d, %d):\t",PX[g0][d],PX[g1][d],PX[g2][d]);
                    for (i = g0; i < g1; i++) {
                        for (j = g1; j < g2+1; j++) {
                            if(PX[g0][d] < PX[g1][d]){
                                DMC[ 0 ][ g ] += DM[ PX[i][d] ][ PX[j][d] ];//+= DM[ cind[i] ][ cind[j] ];
                            }else{
                                DMC[ 0 ][ g ] -= DM[ PX[i][d] ][ PX[j][d] ];//-= DM[ cind[i] ][ cind[j] ];
                            }
                        }
                    }
					//Rprintf("DMC = %d\n",DMC[ 0 ][ g ]);
                }
                
				// do all good steps
				// = apply change to current index and original index IX
				for (g = 0; g < tl[d]; g++) {
                    // do all changes which are an improvement
                    if (DMC[0][g] > 0) {
                        // ll and lr are the numbers of items in the two branches
                        // both are negative iff PX[g1][d] < PX[g0][d]
                        
                        g0 = MAT_ELT(TM,3*g,d,tmp);//TM[ g*3 + 0][d];
                        g1 = MAT_ELT(TM,3*g+1,d,tmp);//TM[ g*3 + 1][d];
                        g2 = MAT_ELT(TM,3*g+2,d,tmp);//TM[ g*3 + 2][d];
                        
                        ll = g1-g0; //PX[g1][d]-PX[g0][d];
                        lr = g2-g1+1; //PX[g2][d]-PX[g1][d]+1;
                        
                        bestcrit += DMC[0][g];
                        if (PX[g1][d] > PX[g0][d]) {
                            for (s = g0; s < g1; s++) {
                                PX[s][d] += lr;
                            }
                            for (s = g1; s < g2+1; s++) {
                                PX[s][d] -= ll;
                            }
                        }else{
                            for (s = g0; s < g1; s++) {
                                PX[s][d] -= lr;
                            }
                            for (s = g1; s < g2+1; s++) {
                                PX[s][d] += ll;
                            }
                        }
                        // Rprintf("node %d:\n",g);
                        //Rprintf("better PX = \n");
                        //for (i = 0; i < tl[d]; i++) {
                        //   Rprintf("%d\t",PX[ i ][d]);
                        //}
                        //anybetter = 1;
                        // an improvement was possible => not yet optimal
                        opt = 0;
                    }
				}
                //if(anybetter == 1){
                
                //rearrange m0 acc. to the new PX order
                // m1 and m0 are still the same => read from m1
                // IX is the old position, PX is the new one ( IX[i][d] == cind[ IX[i][d] ] )
                for (i = 0; i < dimv[d]; i++) {
                    cind[ PX[i][d] ] = IX[i][d];
                }
                //Rprintf("cind = \n");
                //for (i = 0; i < dimv[d]; i++) {
                //	Rprintf("%d\t",cind[ i ]);
                //}
                s1 = cp[d];
                s2 = cpr[d];
                
                for (i = 0; i < dimv[d]; i++) {
                    for (j1 = 0; j1 < s2; j1++) {
                        for (j2 = 0; j2 < s1; j2++) {
                            rind = i*s1 + j1*s1*dimv[d] + j2;
                            //lind = IX[i][d]*s1 + j1*s1*dimv[d] + j2;
                            lind = cind[i]*s1 + j1*s1*dimv[d] + j2;
                            m0[rind] = m1[lind]; 
                        }
                    }
                }
                for (i = 0; i < mln; i++) {
                    m1[i] = m0[i];
                }
                //Rprintf("old PX = \n");
                //for (i = 0; i < dimv[d]; i++) {
                //	Rprintf("%d\t",IX[ i ][d]);
                //}
                //Rprintf("new PX = \n");
                //for (i = 0; i < dimv[d]; i++) {
                //	Rprintf("%d\t",PX[ i ][d]);
                //}
                
				//}
                // update IX (refrain from misusing IX as a temporal PX)
                for (s = 0; s < dimv[d]; s++) {
                    IX[ PX[s][d] ][d] = s;
                }
                
                //Rprintf("new IX = \n");
                //for (i = 0; i < dimv[d]; i++) {
                //    Rprintf("%d\t",IX[ i ][d]);
                //}
                
				imp = best;
                // Rprintf(" m0 = \n");
				//for (i = 0; i < dimv[0]; i++) {
                //   for (j = 0; j < dimv[1]; j++) {
                //     Rprintf("%d\t",m0[i + j*dimv[0]]);
                //}
                //Rprintf("\n\n");
				//}  
			}else{
                // ---------------- free STEP ---------------- //
				best = 0;
		//		Rprintf("GO!\n\n");
				for (i = 0; i < dimv[d]-1; i++) {
					for (j = i+1; j < dimv[d]; j++) {
						DMCl[ i ][ j ] = DMCl[ i ][ j-1 ] - (j-i)*DM[ cind[i] ][ cind[j] ]; //+ DM[ cind[i] ][ cind[j] ];
                    }
                }
                
                for (i = 0; i < dimv[d]-1; i++) {
					for (j = i+1; j < dimv[d]; j++) {
						DMCw[ i ][ j ] = DMCw[ i ][ j-1 ] + DM[ cind[j] ][ cind[i] ]; //+ DM[ cind[i] ][ cind[j] ];
					}
                    for (j = i+1; j < dimv[d]; j++) {
						DMCw[ i ][ j ] += DMCw[ i ][ j-1 ]; //+ DM[ cind[i] ][ cind[j] ];
                       
					}
                   
				}
                
                for (i = 0; i < dimv[d]-1; i++) {
					for (j = i+1; j < dimv[d]; j++) {
						DMC[ i ][ j ] = DMCl[ i ][ j ] + DMCw[ i ][ j ]; 
                        if(DMC[ i ][ j ] > best){
                            ii = i;
                            jj = j;
                            best = DMC[ i ][ j ];
                        }
                        
					}
                    
				}
               
                
                
                for (i = 1; i < dimv[d]; i++) {
					for (j = i-1; j >= 0; j--) {
						DMCl[ i ][ j ] = DMCl[ i ][ j+1 ] - (i-j)*DM[ cind[j] ][ cind[i] ]; //+ DM[ cind[i] ][ cind[j] ];
					}
				}
                for (i = 1; i < dimv[d]; i++) {
					for (j = i-1; j >= 0; j--) {
						DMCw[ i ][ j ] = DMCw[ i ][ j+1 ] + DM[ cind[i] ][ cind[j] ]; //+ DM[ cind[i] ][ cind[j] ];
					}
                    for (j = i-1; j >= 0; j--) {
						DMCw[ i ][ j ] += DMCw[ i ][ j+1 ]; //+ DM[ cind[i] ][ cind[j] ];
					}
				}
              
				for (i = 1; i < dimv[d]; i++) {
					for (j = i-1; j >= 0; j--) {
						DMC[ i ][ j ] = DMCl[ i ][ j ] + DMCw[ i ][ j ]; 
					}
				}
                
                for (s = 0; s < dimv[d]-1; s++) {
                    for (t=s+1; t < dimv[d]; t++) {
                        // i moves between s and t from the left
                        if (s > 0) {
                        for (i = 0; i < s; i++) {
                            for (j = s; j < t; j++) {
                                DMC[i][j] = DMC[i][j] + DM[ cind[s] ][ cind[t] ] - DM[ cind[t] ][ cind[s] ];
                            }
                        }
                        }
                        // i moves between s and t from the right
                        if( t < dimv[d]-1 ){
                        for (i = dimv[d]-1; i > t; i--) {
                            for (j = t; j > s; j--) {
                                DMC[i][j] = DMC[i][j] + DM[ cind[s] ][ cind[t] ] - DM[ cind[t] ][ cind[s] ];
                            }
                        }
                        }
                        // i moves from in between s and t out to the left of s or the right if t
                        if(s+1 < t){
                        for (i = s+1; i < t; i++) {
                            for (j = 0; j <= s; j++) {
                                DMC[i][j] = DMC[i][j] - DM[ cind[s] ][ cind[t] ] + DM[ cind[t] ][ cind[s] ];
                            }
                            for (j = t; j < dimv[d]; j++) {
                                DMC[i][j] = DMC[i][j] - DM[ cind[s] ][ cind[t] ] + DM[ cind[t] ][ cind[s] ];
                            }
                        }
                        }
						// t = i
						// <- i
						if(s+1<t){
						for (j = s+1; j < t; j++) {
							DMC[t][j] = DMC[t][j] - (t-j)*DM[ cind[s] ][ cind[t] ] + (t-j)*DM[ cind[t] ][ cind[s] ];
						}
						}
						
						// i ->
						if(t+1 < dimv[d]){
						for (j = t+1; j < dimv[d]; j++) {
							DMC[t][j] = DMC[t][j] + (j-t)*DM[ cind[s] ][ cind[t] ] - (j-t)*DM[ cind[t] ][ cind[s] ];
						}
						}
						// s = i
						// i ->
						if(s+1 < t){
						for (j = s+1; j < t; j++) {
							DMC[s][j] = DMC[s][j] - (t-j)*DM[ cind[s] ][ cind[t] ] + (t-j)*DM[ cind[t] ][ cind[s] ];
						}
						}
						
						// <- i
						if(s > 0){
						for (j = 0; j < s; j++) {
							DMC[s][j] = DMC[s][j] + (s-j)*DM[ cind[s] ][ cind[t] ] - (s-j)*DM[ cind[t] ][ cind[s] ];
						}
						}
                    }
                }
                for (i = 1; i < dimv[d]; i++) {
					for (j = i-1; j >= 0; j--) {
						if(DMC[ i ][ j ] > best){
                            ii = i;
                            jj = j;
                            best = DMC[ i ][ j ];
                        }
					}
				}
                
				//Rprintf("DMCw\n");
				//for (i = 0; i < dimv[d]; i++) {
				//	for (j = 0; j < dimv[d]; j++) {
				///		Rprintf("%d\t",DMCw[ i ][ j ]);
				//	}
                //    Rprintf("\n");
				//}
               // Rprintf("\n");
				//Rprintf("DMCl\n");
				//for (i = 0; i < dimv[d]; i++) {
				//	for (j = 0; j < dimv[d]; j++) {
				////	}
                 ///   Rprintf("\n");
				//}
		//Rprintf("\n");
  //  Rprintf("DMC\n");
				//for (i = 0; i < dimv[d]; i++) {
					//for (j = 0; j < dimv[d]; j++) {
					//	Rprintf("%d\t",DMC[ i ][ j ]);
					//}
                   // Rprintf("\n");
				//}
               // Rprintf("\n");
				
				
				
				// do best step
				// = apply change to current index and original index IX
				if(best > eps){
					// Rprintf(" ---- ---- %d / %d is better by %f\t", ii,jj, best);
					bestcrit += best;
					tmp = IX[ii][d];
					tmp2 = cind[ii];
					if (ii < jj) {
						for (s = ii; s < jj; s++) {
							IX[s][d] = IX[s+1][d];
							cind[s] = cind[s+1];
						}
						IX[jj][d] = tmp;
						cind[jj] = tmp2;
					}
					if (ii > jj) {
						for (s = ii; s > jj; s--) {
							IX[s][d] = IX[s-1][d];
							cind[s] = cind[s-1];
						}
						IX[jj][d] = tmp;
						cind[jj] = tmp2;
					}
					// an improvement was possible => not yet optimal
					opt = 0;
				}else{
					// Rprintf("### -------- stop ------ #\n");
					
					//rearrange m0 acc. to the new IX order
					// m1 and m0 are still the same => read from m1
					s1 = cp[d];
					s2 = cpr[d];
					
					for (i = 0; i < dimv[d]; i++) {
						for (j1 = 0; j1 < s2; j1++) {
							for (j2 = 0; j2 < s1; j2++) {
								rind = i*s1 + j1*s1*dimv[d] + j2;
								//lind = IX[i][d]*s1 + j1*s1*dimv[d] + j2;
								lind = cind[i]*s1 + j1*s1*dimv[d] + j2;
								m0[rind] = m1[lind]; 
							}
						}
					}
					for (i = 0; i < mln; i++) {
						m1[i] = m0[i];
					}
				}
				Rprintf("new IX = \n");
                for (i = 0; i < dimv[d]; i++) {
                    Rprintf("%d\t",IX[ i ][d]);
                }
				Rprintf("free %d\n\n", d);
				
				
				imp = best;
				Rprintf("imp %f\n\n", imp);
				
				// Rprintf(" m0 = \n");
				//for (i = 0; i < dimv[0]; i++) {
				//    for (j = 0; j < dimv[1]; j++) {
				//      Rprintf("%d\t\t",m0[i + j*dimv[0]]);
                //}
                //Rprintf("\n");
				//}
				
			}
			
		}
		/// ''''' END INTRA STEP 2 ''''' ///
		
        /// -----------------------------------------------------------------------				
		// decide which d is next
		d = (d+1)%nd;
		
		//
        //  Rprintf("\n new order %d\n",d);
		newcrit = mvclasscrit2(mv, IXp, dimv, nd);
        //Rprintf("newcrit = %f\n",newcrit);
		//Rprintf("bestcrit = %f\n",bestcrit);
        //if(newcrit <= bestcrit){
        //   opt = 1;
        //}else{
        //    bestcrit = newcrit;
        // }
		
	}
	/// ''''' END MAIN STEP ''''' ///
	free(TM);
	
	SEXP out = allocVector(REALSXP, biglen + 1);
	
	k = 0;
	for(i = 0; i < nd; i++){
		for (j = 0; j < dimv[i]; j++) {
			REAL(out)[k] = IX[j][i];
			k++;
		}
	}
	REAL(out)[biglen] = bestcrit;
	Rprintf("bestcrit = %f\n",bestcrit);
	return out;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //




// --------------------------------------------------------------------------------------- //
// -----------------------------         QUICK FECHNER         --------------------------- //
// --------------------------------------------------------------------------------------- //

SEXP quickfechner(SEXP M, SEXP dims, SEXP vs, SEXP exz){
    
	
	int ex = INTEGER(exz)[0];
	int i,j,k;//, d, dr, dc;
	float  dvr, dvc;//dv
	//int nd = LENGTH(dims);
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
    
	float X[n][m];// dissimilarity matrix
	short C[n][m];// change indicator matrix
	short P[n][m];// parent matrix
	//short C2[n][m];
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//		X[i][j] = REAL(M)[n*j + i];
	//		C[i][j] = 0;
	//		P[i][j] = i;
	
	//	}
	//}
	
	if(ex == 1){
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				X[i][j] = REAL(M)[n*j + i];
				if(X[i][j] == 0){
					C[i][j] = 0;
				}else{
					C[i][j] = 1;
				}
				P[i][j] = i;
			}
		}
	}else{
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				X[i][j] = REAL(M)[n*j + i];
				C[i][j] = 1;
				P[i][j] = i;
			}
		}
	}
	
	int optimal = 0;

	if (INTEGER(vs)[0] == 0) { // additive algorithm
	//initial step
		//for (i = 0; i < n; i++) {
			//for (j = 0; j < m; j++) {
			//	for (k = 0; k < m; k++) {
					
					
				//	// do not consider zero-entries (i.e. nodes are not connected)
				//	if( X[i][k]*X[k][j] == 0 & ex == 1){
				//		dv = X[i][j]+1;
				//	}else{
				//		dv = X[i][k] +  X[k][j];
				//	}
					
				//	if (X[i][j] > dv ) {
				//		//Rprintf("%d, %d ,%d\n",i,j,k);
				//		X[i][j] = dv;
				//		C[i][j] = 1;
				//		P[i][j] = k;
				//	}
				//}
			//}
		//}

		//Rprintf("\n start loop \n\n");
	while(optimal ==  0){
		optimal = 1;
	//further steps
	for (i = 0; i < n; i++) {
		for (k = 0; k < m; k++) {
			if(C[i][k] == 1){ // [i][j] ??
				for (j = 0; j < m; j++) {
					// i j is either X[i][k] or X[k][j] from above
					
					
					
					// do not consider zero-entries (i.e. nodes are not connected)
					if( (X[i][k]*X[j][i] == 0) & (ex == 1)){
						dvr = X[j][k]+1;
					}else{
						dvr = X[i][k] +  X[j][i];
					}	
					if( (X[i][k]*X[k][j] == 0) & (ex == 1)){
						dvc = X[i][j]+1;
					}else{
						dvc = X[i][k] +  X[k][j];
					}
					
					if (X[j][k] > dvr ) {
						X[j][k] = dvr;
						C[j][k] = 1;
						P[j][k] = P[i][k];
						optimal = 0;
					}
					if (X[i][j] > dvc ) {
						X[i][j] = dvc;
						C[i][j] = 1;
						P[i][j] = P[k][j];//k;
						optimal = 0;
					}
				}
				C[i][k] = 0;
			}
		}
	}
		
	}
    }else{ // multiplicative algorithm
        //initial step
		//for (i = 0; i < n; i++) {
			//for (j = 0; j < m; j++) {
			//	for (k = 0; k < m; k++) {
					
			//		// do not consider zero-entries (i.e. nodes are not connected)
				//	dv = X[i][k] *  X[k][j];
					
				//	if( dv == 0 & ex == 1){
				//		dv = X[i][j]-1;
				//	}
					
				//	if (X[i][j] < dv ) {
				//		//Rprintf("%d, %d ,%d\n",i,j,k);
				//		X[i][j] = dv;
				//		C[i][j] = 1;
				//		P[i][j] = k;
				//	}
				//}
			//}
		//}
        
		//Rprintf("\n start loop \n\n");
        while(optimal ==  0){
            optimal = 1;
            //further steps
            for (i = 0; i < n; i++) {
                for (k = 0; k < m; k++) {
                    if(C[i][k] == 1){
                        for (j = 0; j < m; j++) {
                            // i j is either X[i][k] or X[k][j] from above
                          
							// do not consider zero-entries (i.e. nodes are not connected)
							dvr = X[i][k] *  X[j][i];
							dvc = X[i][k] *  X[k][j];
							
							if( (dvr == 0) & (ex == 1)){
								dvr = X[j][k]-1;
							}	
							if( (dvc == 0) & (ex == 1)){
								dvc = X[i][j]-1;
							}
							
                            if (  X[j][k] < dvr ) {
                                X[j][k] = dvr;
                                C[j][k] = 1;
								P[j][k] = P[i][k];
                                optimal = 0;
                            }
                            if (X[i][j] < dvc ) {
                                X[i][j] = dvc;
                                C[i][j] = 1;
								P[i][j] = P[k][j];//k;
                                optimal = 0;
                            }
                        }
                        C[i][k] = 0;
                    }
                }
            }
            
        }
    
    }
	SEXP out = allocVector(REALSXP, 2*n*m);
	k = 0;
	for (j = 0; j < m; j++) {
	for(i = 0; i < n; i++){
		
			//Rprintf("%f\t",X[i][j]);
			REAL(out)[k] = X[i][j];
			k++;
		}
		//Rprintf("\n");
	}
	
	for (j = 0; j < m; j++) {
	for(i = 0; i < n; i++){
		
			//Rprintf("%f\t",X[i][j]);
			REAL(out)[k] = P[i][j]+1;
			k++;
		}
		//Rprintf("\n");
	}
	
	return out;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //







// --------------------------------------------------------------------------------------- //
// -----------------------------         QUICK FECHNER  II       --------------------------- //
// --------------------------------------------------------------------------------------- //
SEXP quickfechner2(SEXP M, SEXP dims, SEXP vs, SEXP exz){
    
	
	int ex = INTEGER(exz)[0];
	int i,j,k, d;// dr, dc;
	float  dvr, dvc;//dv
	//int nd = LENGTH(dims);
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
    
	float X[n][m];// dissimilarity matrix
	short C[n][m];// change indicator matrix
	short P[n][m];// parent matrix
	//short C2[n][m];
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//		X[i][j] = REAL(M)[n*j + i];
	//		C[i][j] = 0;
	//		P[i][j] = i;
	
	//	}
	//}
	
	if(ex == 1){
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				X[i][j] = REAL(M)[n*j + i];
				if(X[i][j] == 0){
					C[i][j] = 0;
				}else{
					C[i][j] = 1;
				}
				P[i][j] = i;
			}
		}
	}else{
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				X[i][j] = REAL(M)[n*j + i];
				C[i][j] = 1;
				P[i][j] = i;
			}
		}
	}
	
	int optimal = 0;
	
	if (INTEGER(vs)[0] == 0) { // additive algorithm
		//initial step
		//for (i = 0; i < n; i++) {
		//for (j = 0; j < m; j++) {
		//	for (k = 0; k < m; k++) {
		
		
		//	// do not consider zero-entries (i.e. nodes are not connected)
		//	if( X[i][k]*X[k][j] == 0 & ex == 1){
		//		dv = X[i][j]+1;
		//	}else{
		//		dv = X[i][k] +  X[k][j];
		//	}
		
		//	if (X[i][j] > dv ) {
		//		//Rprintf("%d, %d ,%d\n",i,j,k);
		//		X[i][j] = dv;
		//		C[i][j] = 1;
		//		P[i][j] = k;
		//	}
		//}
		//}
		//}
		
		//Rprintf("\n start loop \n\n");
		while(optimal ==  0){
			optimal = 1;
			//further steps
			for(d = 1; d < n-1; d++){
			for (i = 0; i < n; i++) {
				for (k = 0; k < m; k++) {
					if(C[i][k] == 1){ // [i][j] ??
						for (j = mymax2(i-d,0); j < mymin2(i+d,n); j++) {
							// i j is either X[i][k] or X[k][j] from above
							// do not consider zero-entries (i.e. nodes are not connected)
								
							if( (X[i][k]*X[k][j] == 0) & (ex == 1)){
								dvc = X[i][j]+1;
							}else{
								dvc = X[i][k] +  X[k][j];
							}
							if (X[i][j] > dvc ) {
								X[i][j] = dvc;
								C[i][j] = 1;
								P[i][j] = P[k][j];//k;
								optimal = 0;
							}
						}
						for (j = mymax2(k-d,0); j < mymin2(k+d,n); j++) {	
							if( (X[i][k]*X[j][i] == 0) & (ex == 1)){
								dvr = X[j][k]+1;
							}else{
								dvr = X[i][k] +  X[j][i];
							}
							if (X[j][k] > dvr ) {
								X[j][k] = dvr;
								C[j][k] = 1;
								P[j][k] = P[i][k];
								optimal = 0;
							}
						}
						C[i][k] = 0;
					}
				}
			}}
			
		}
    }else{ // multiplicative algorithm
        //initial step
		//for (i = 0; i < n; i++) {
		//for (j = 0; j < m; j++) {
		//	for (k = 0; k < m; k++) {
		
		//		// do not consider zero-entries (i.e. nodes are not connected)
		//	dv = X[i][k] *  X[k][j];
		
		//	if( dv == 0 & ex == 1){
		//		dv = X[i][j]-1;
		//	}
		
		//	if (X[i][j] < dv ) {
		//		//Rprintf("%d, %d ,%d\n",i,j,k);
		//		X[i][j] = dv;
		//		C[i][j] = 1;
		//		P[i][j] = k;
		//	}
		//}
		//}
		//}
        
		//Rprintf("\n start loop \n\n");
        while(optimal ==  0){
            optimal = 1;
            //further steps
            for (i = 0; i < n; i++) {
                for (k = 0; k < m; k++) {
                    if(C[i][k] == 1){
                        for (j = 0; j < m; j++) {
                            // i j is either X[i][k] or X[k][j] from above
							
							// do not consider zero-entries (i.e. nodes are not connected)
							dvr = X[i][k] *  X[j][i];
							dvc = X[i][k] *  X[k][j];
							
							if( (dvr == 0) & (ex == 1)){
								dvr = X[j][k]-1;
							}	
							if( (dvc == 0) & (ex == 1)){
								dvc = X[i][j]-1;
							}
							
                            if (  X[j][k] < dvr ) {
                                X[j][k] = dvr;
                                C[j][k] = 1;
								P[j][k] = P[i][k];
                                optimal = 0;
                            }
                            if (X[i][j] < dvc ) {
                                X[i][j] = dvc;
                                C[i][j] = 1;
								P[i][j] = P[k][j];//k;
                                optimal = 0;
                            }
                        }
                        C[i][k] = 0;
                    }
                }
            }
            
        }
		
    }
	SEXP out = allocVector(REALSXP, 2*n*m);
	k= 0;
	for (j = 0; j < m; j++) {
		for(i = 0; i < n; i++){
			
			//Rprintf("%f\t",X[i][j]);
			REAL(out)[k] = X[i][j];
			k++;
		}
		//Rprintf("\n");
	}
	
	for (j = 0; j < m; j++) {
		for(i = 0; i < n; i++){
			
			//Rprintf("%f\t",X[i][j]);
			REAL(out)[k] = P[i][j]+1;
			k++;
		}
		//Rprintf("\n");
	}
	
	return out;
}





// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //







float indcrit2(int *m0, int *rord, int *cord, int n, int m){
    // hamming distance criterion	
	
	int i,j;
	//int (*MX)[m] = (int (*)[m])m0;
    
	
	float loss = 0.0;
	
	//float MY[n][m];
	long M2[n][m]; // cumsum by column
	long M1[n][m]; // cumsum by row
	
	//set the first row and the last column(s)
	for( i = 0; i < n; i++ ){ 
		M1[i][m-2] = m0[rord[i] + cord[m-1]*n];
		M1[i][m-1] = 0;
	}
	for( j = 0; j < m; j++ ){ 
		M2[0][j] = 0;
		M2[1][j] = m0[rord[0] + cord[j]*n];
	}
	
	
	// compute M1
	if(m > 2){ 
		for( j = m - 3; j >= 0 ; j-- ){
			for( i=0; i < n; i++ ){
				M1[i][j] = m0[rord[i] + cord[j+1]*n] + M1[i][j+1];
			}
		}
	}
	for( j = m - 3; j >= 0 ; j-- ){
		for( i = 0; i < n; i++ ){
			M1[i][j] = M1[i][j] + M1[i][j+1];
		}
	}
	
	for( i = 1; i < n; i++ ){
		for( j= m - 2; j >= 0 ; j-- ){
			M1[i][j] = M1[i][j] + M1[i-1][j];
		}
	}
	
	// compute M2	
	if(n > 2){
		for( i=2; i < n; i++ ){
			for( j = 0; j < m; j++ ){
				M2[i][j] = m0[rord[i-1] + cord[j]*n] + M2[i-1][j];
			}
		}
	}
	for( i=2; i < n; i++ ){
		for( j = 0; j < m; j++ ){
			M2[i][j] = M2[i][j] + M2[i-1][j];
		}
	}
	
	for( j = m - 2; j >= 0 ; j-- ){
		for( i = 1; i < n; i++ ){
			M2[i][j] = M2[i][j] + M2[i][j+1];
		}
	}
    	// compute the criterion
	for( i = 0; i < n; i++ ){
		for( j = 0; j < m; j++ ){
			loss = loss + (M1[i][j] + M2[i][j]) * m0[rord[i] + cord[j]*n];
		}
	}
	return loss;
}



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //



// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //



SEXP quickhamm2d(SEXP M, SEXP dims, SEXP pv){
	
	// M is byrow
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int px = INTEGER(pv)[1];
	int py = INTEGER(pv)[0];
	
	
	
	//Rprintf("%d\t%d\t%d\t%d\t%d\t",n,m,px,py);
	
	int i,j,s,tmp, bi, bj;
	
	
	// index vectors with initial values 0..n-1 and 0..m-1 containing the actual orders
	int RI[n];
	int CI[m];
    int TRI[n];
	int TCI[m];
	for (i = 0; i < n; i++) {
		RI[i] = i;
        TRI[i] = i;
	}
	for (j = 0; j < m; j++) {
		CI[j] = j;
        TCI[j] = j;
	}
	
	// Delta matrices for rows and columns
	// crit( i -> i2 ) - crit( i2 -> i ) if i < i2 and -(...) else
	// i.e. this is 'lost' when putting i to i2
	
	
	
	
    
	
	//Rprintf("check02");
	// the criterion
	
    float currcrit = 0;
	
	
	//an marray pointer to M
	
    int* mv = &INTEGER(M)[0];
     int* pri = &TRI[0];  
    int* pci = &TCI[0]; 
    
	float bestcrit = indcrit2(mv, pri,pci,n,m);
    int opt = 0;
    while (opt == 0) {
        opt = 1;
        if(px == 1){
        bi = n+m;
        bj = n+m;
        //try column i to pos j
        for (i=0; i < m; i++) {
            for (j = 0; j < m; j++) {
                if( j > i ){
                    tmp = TCI[i];
                    for (s = i; s < j; s++) {
                        TCI[s] = TCI[s+1];
                    }
                    TCI[j] = tmp;
                }
                if( i > j ){
                    tmp = TCI[i];
                    for (s = i; s > j; s--) {
                        TCI[s] = TCI[s-1];
                    }
                    TCI[j] = tmp;
                }
                currcrit = indcrit2(mv, pri, pci, n, m );
                if(currcrit < bestcrit){
                    bestcrit = currcrit;
                    bi = i;
                    bj = j;
                }
                for (s = 0; s < m; s++) {
                    TCI[s] = CI[s];
                }
            }
        }
        if(bi < n+m){
             opt = 0;
            j = bj;
            i = bi;
            if( j > i ){
                tmp = CI[i];
                for (s = i; s < j; s++) {
                    CI[s] = CI[s+1];
                }
                CI[j] = tmp;
            }
            if( i > j ){
                tmp = CI[i];
                for (s = i; s > j; s--) {
                    CI[s] = CI[s-1];
                }
                CI[j] = tmp;
            }
            for (s = 0; s < m; s++) {
                TCI[s] = CI[s];
            }
        }
        }
        if(py == 1){
            bi = n+m;
            bj = n+m;
        //try row i to pos j
        for (i=0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if( j > i ){
                    tmp = TRI[i];
                    for (s = i; s < j; s++) {
                        TRI[s] = TRI[s+1];
                    }
                    TRI[j] = tmp;
                }
                if( i > j ){
                    tmp = TRI[i];
                    for (s = i; s > j; s--) {
                        TRI[s] = TRI[s-1];
                    }
                    TRI[j] = tmp;
                }
                currcrit = indcrit2(mv, pri, pci, n, m );
                if(currcrit < bestcrit){
                    bestcrit = currcrit;
                    bi = i;
                    bj = j;
                    opt = 0;
                }
                for (s = 0; s < n; s++) {
                    TRI[s] = RI[s];
                }
            }
        }
        if(bi < n+m){
            opt = 0;
            j = bj;
            i = bi;
            if( j > i ){
                tmp = RI[i];
                for (s = i; s < j; s++) {
                    RI[s] = RI[s+1];
                }
                RI[j] = tmp;
            }
            if( i > j ){
                tmp = RI[i];
                for (s = i; s > j; s--) {
                    RI[s] = RI[s-1];
                }
                RI[j] = tmp;
            }
            bi = n+m;
            bj = n+m;
            for (s = 0; s < n; s++) {
                TRI[s] = RI[s];
            }
        }
        }
    }
	
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		REAL(out)[i] = RI[i]; //+1
	}
	for( j = 0; j < m;j++ ){
		REAL(out)[n+j] = CI[j]; //+1
	}
	
	REAL(out)[n+m] = bestcrit;
	//Rprintf("\nOptimization ended after %d steps with crtv=%f\n",steps,crtv);
	return out;
}


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //



// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //

SEXP quicktileMR(SEXP M, SEXP dims, SEXP pv, SEXP brks, SEXP iter){
	
	// M is byrow
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int px = INTEGER(pv)[1];
	int py = INTEGER(pv)[0];
	
	int nb = LENGTH(brks);
	
	
	
	//Rprintf("%d\t%d\t%d\t%d\t%d\t",n,m,px,py);
	
	int i,j,i2,j2,k,tmp;
	
	
	// index vectors with initial values 0..n-1 and 0..m-1 containing the actual orders
	int RI[n];
	int CI[m];
	for (i = 0; i < n; i++) {
		RI[i] = i;
	}
	for (j = 0; j < m; j++) {
		CI[j] = j;
	}
	int idv[m];
	
	int ii,ii2;//,jj,jj2;
	
	// Delta matrices for rows and columns
	// crit( i -> i2 ) - crit( i2 -> i ) if i < i2 and -(...) else
	// i.e. this is 'lost' when putting i to i2
	
	
   
	// Sum of delta matrices 
	

	float *DC2 = calloc(m*m,sizeof(float));
	
	float *DR2 = calloc(n*n,sizeof(float));
	float *DC = calloc(m*m,sizeof(float));
	
	float *DR = calloc(n*n,sizeof(float));
	
	int *MS = calloc(n*m,sizeof(int));	
	int *MT = calloc(n*m,sizeof(int));	
	
	
	
	if (py == 1) {
		for (i = 0; i < n; i++) {
			MAT_ELT(DR2,i,i,n) = 0;
			for (i2 = 0; i2 < n; i2++) {
				MAT_ELT(DR,i,i2,n) = 0;
			}
		}
	}
	if (px == 1) {
		for (j = 0; j < m; j++) {
			MAT_ELT(DC2,j,j,m) = 0;
			
			for (j2 = 0; j2 < m; j2++) {
				MAT_ELT(DC,j,j2,m) = 0;
			}
		}
	}
	
	//Rprintf("check02");
	// the criterion
	float crtv = 0;
	float td;
	
	//an marray pointer to M
	
	//int* mv = &INTEGER(M)[0];
	//int (*MT)[m] = (int (*)[m])mv;
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//	Rprintf("%d\t", MT[i][j]);
	//}
	//Rprintf("\n");
	//}
	
	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++) {
			MAT_ELT(MT,i,j,n) = INTEGER(M)[i+j*n];
		}
	}
	int ding = 0;
	int ik = 0;
	
	
	float bestc = 0;
	float bestr = 0;
	float best = 0;
	int r1 = 0;
	int r2 = 0;
	int c1 = 0;
	int c2 = 0;
	int opt = 0;
	int steps = 0;
	int* id;
	int bl;
	
	while (ik < INTEGER(iter)[0]) {
		ik++;
	for (bl = 0+ding; bl < nb+ding; bl++) {
		id = &idv[0];
	
	
	// --------------------- compute the initial delta values --------------------- //
	
	// cumsum by columns
	for (j = 0; j < m; j++) {
		MAT_ELT(MS,0,j,n) = MAT_ELT(MT,0,j,n);
	}
	for (i = 1; i < n; i++) {
		for (j = 0; j < m; j++) {
			MAT_ELT(MS,i,j,n) = MAT_ELT(MT,i,j,n) + MAT_ELT(MS,i-1,j,n);
		}
	}
	
	// delta C
	
	for (j = 0; j < m-1; j++) {
		for (j2 = j+1; j2 < m; j2++) {
			for (i = 1; i < n; i++) {
				td = MAT_ELT(MT,i,j,n) * MAT_ELT(MS,i-1,j2,n);
				crtv += td;
				//DC[j][j2] += ( MT[i][j]*MS[i-1][j2] - MT[i][j2]*MS[i-1][j] ); // loss
				MAT_ELT(DC,j,j2,m) += ( td - MAT_ELT(MT,i,j2,n) * MAT_ELT(MS,i-1,j,n) ); // loss
			}
			MAT_ELT(DC,j2,j,m) = - MAT_ELT(DC,j,j2,m);
		}
	}
	// also computes the initial criterion
	
	if (py == 1) {
		// cumsum by rows
		for (i = 0; i < n; i++) {
			MAT_ELT(MS,i,0,n) = MAT_ELT(MT,i,0,n);
		}
		for (i = 0; i < n; i++) {
			for (j = 1; j < m; j++) {
				MAT_ELT(MS,i,j,n) = MAT_ELT(MT,i,j,n) + MAT_ELT(MS,i,j-1,n);
			}
		}
		
		// delta R
		for (i = 0; i < n-1; i++) {
			for (i2 = i+1; i2 < n; i2++) {
				for (j = 1; j < m; j++) {
					MAT_ELT(DR,i,i2,n) += ( MAT_ELT(MT,i,j,n) * MAT_ELT(MS,i2,j-1,n) - MAT_ELT(MT,i2,j,n) * MAT_ELT(MS,i,j-1,n) ); // loss
				}
				MAT_ELT(DR,i2,i,n) = - MAT_ELT(DR,i,i2,n);
			}
		}
	}
	
	
	//Rprintf("\n\nstart crtv = %f\n\n",crtv);
	
	// --------------------- start the main algorithm --------------------- //	
	
	
	// look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to i2)
	// 'best' is the decrease in crtv, r1, r2, c1, c2 are the corresp. indices

		bestc = 0;
		bestr = 0;
		best = 0;
		r1 = 0;
		r2 = 0;
		c1 = 0;
		c2 = 0;
		opt = 0;
		steps = 0;
	
	while( opt == 0 ){	
		
		
		//look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to k)
		
		if (py == 1) {
			for (i = 0; i < n-1; i++) {
				ii = id[i];
				
				for (i2 = i+1; i2 < n; i2++) {
					ii2 = id[i2];
					
					MAT_ELT(DR2,i,i2,n) = MAT_ELT(DR2,i,i2-1,n)+ MAT_ELT(DR, RI[i] , RI[i2] ,n);
					if (MAT_ELT(DR2,i,i2,n) > bestr) {
						
						bestr = MAT_ELT(DR2,i,i2,n);
						r1 = i;
						r2 = i2;
					}
				}
				
			}
			
			for (i = 1; i < n; i++) {
				for (i2 = i-1; i2 >= 0; i2--) {
					MAT_ELT(DR2,i,i2,n) = MAT_ELT(DR2,i,i2+1,n)+  MAT_ELT(DR, RI[i2] , RI[i] ,n);
					if (MAT_ELT(DR2,i,i2,n) > bestr) {
						
						bestr = MAT_ELT(DR2,i,i2,n);
						r1 = i;
						r2 = i2;
					}
				}
				
			}
		}
		
		
		if (px == 1) {
			for (j = 0; j < m-1; j++) {
				for (j2 = j+1; j2 < m; j2++) {
					MAT_ELT(DC2,j,j2,m) = MAT_ELT(DC2,j,j2-1,m)+ MAT_ELT(DC, CI[j] , CI[j2] ,m);
					if (MAT_ELT(DC2,j,j2,m) > bestc) {
						
						bestc = MAT_ELT(DC2,j,j2,m);
						c1 = j;
						c2 = j2;
					}
				}
				
			}
			
			for (j = 1; j < m; j++) {
				for (j2 = j-1; j2 >= 0; j2--) {
					MAT_ELT(DC2,j,j2,m) = MAT_ELT(DC2,j,j2+1,m)+ MAT_ELT(DC, CI[j2] , CI[j] ,m);
					if (MAT_ELT(DC2,j,j2,m) > bestc) {
						
						bestc = MAT_ELT(DC2,j,j2,m);
						c1 = j;
						c2 = j2;
					}
				}
				
			}
		}
		
		
		
		// two columns OR two rows are exchanged. The DR/DC values will change whilst DC/DR remains unchanged.
		// the original data matrix and the delta matrices DC and DR and MS remain unchanged with repect to their row/column orders.
		// BUT the values change dependent on the permutations! Hence the sums (DR2 and DC2) use CI[j] and RI[i]
		
		//Rprintf("\n\n r1 = %d, r2 = %d, bestr = %f\n",r1,r2, bestr);
		//Rprintf("\n\n c1 = %d, c2 = %d, bestc = %f\n",c1,c2, bestc);
		
		//bestc = bestr+1;
		
		
		// reset bestc and bestr
		if( bestc > 0 || bestr > 0 ){	
			
			if (bestc > bestr) {
				best = bestc;
				
				if (c1 < c2) {
					// cumulative partial rowsums
					for (i = 0; i < n; i++) {
						MAT_ELT(MS,i, CI[c2] ,n) = MAT_ELT(MT,i, CI[c2] ,n);
					}
					if (c1+1 < c2) { // not neighboring
						for (k = c2-1; k > c1; k--) {
							for (i = 0; i < n; i++) {
								MAT_ELT(MS,i, CI[k] ,n) = MAT_ELT(MT,i, CI[k] ,n) + MAT_ELT(MS,i, CI[k+1] ,n);
							}
						}
					}
					if (py == 1) {
						// changes to DR
						for (i = 0; i < n-1; i++) {
							for (i2 = i+1; i2 < n; i2++) {
								//							Rprintf("c1 = %d, CI[c1+1] = %d\n",c1,CI[c1+1]);
								//							Rprintf("----------\n");
								//							Rprintf("MT[%d,%d] = %d, MS[%d,%d] = %d, prod = %f\n",i2,CI[c1],MT[i2][ CI[c1] ],i,CI[c1+1],MS[i][ CI[c1+1] ],2*(float)MT[i2][ CI[c1] ]*MS[i][ CI[c1+1] ]);
								//							Rprintf("----------\n");
								//							Rprintf("\nMT[%d,%d] = %d, MS[%d,%d] = %d, prod = %f\n",i,CI[c1],MT[i][ CI[c1] ],i2,CI[c1+1],MS[i2][ CI[c1+1] ],2*(float)MT[i][ CI[c1] ]*MS[i2][ CI[c1+1] ]);
								//Rprintf("----------\n");
								MAT_ELT(DR,i,i2,n) -= 2*(  (float)MAT_ELT(MT,i2, CI[c1] ,n) * MAT_ELT(MS,i, CI[c1+1] ,n) -  (float)MAT_ELT(MT,i, CI[c1] ,n) * MAT_ELT(MS,i2, CI[c1+1] ,n) );
								MAT_ELT(DR,i2,i,n) = -MAT_ELT(DR,i,i2,n);
							}
						}
					}
					
					// index vector CI changes
					tmp = CI[c1];
					for (i = c1; i < c2; i++) {
						CI[i] = CI[i+1];
					}
					CI[c2] = tmp;
					//Rprintf("------------M-------------\n");
					
					//for (i = 0; i < n; i++) {
					//	for (j = 0; j < m; j++) {
					//		Rprintf("%d\t",MT[ RI[i] ][ CI[j] ]);
					//	}
					//	Rprintf("\n");
					//}
					
					//Rprintf("----------------------------\n");
					
				}else {
					// cumulative partial rowsums
					for (i = 0; i < n; i++) {
						MAT_ELT(MS,i, CI[c2] ,n) = MAT_ELT(MT,i, CI[c2] ,n);
					}
					if (c2+1 < c1) { // not neighboring
						for (k = c2+1; k < c1; k++) {
							for (i = 0; i < n; i++) {
								MAT_ELT(MS,i, CI[k] ,n) = MAT_ELT(MT,i, CI[k] ,n) + MAT_ELT(MS,i, CI[k-1] ,n);
							}
						}
					}
					if (py == 1) {
						// changes to DR
						for (i = 0; i < n-1; i++) {
							for (i2 = i+1; i2 < n; i2++) {
								MAT_ELT(DR,i,i2,n) -= 2*(  (float)MAT_ELT(MT,i, CI[c1] ,n) * MAT_ELT(MS,i2, CI[c1-1] ,n) - (float)MAT_ELT(MT,i2, CI[c1] ,n) * MAT_ELT(MS,i, CI[c1-1] ,n) ); // neg. of c1 < c2
								MAT_ELT(DR,i2,i,n) = -MAT_ELT(DR,i,i2,n);
							}
						}
					}
					
					// index vector CI changes
					tmp = CI[c1];
					for (i = c1; i > c2; i--) {
						CI[i] = CI[i-1];
					}
					CI[c2] = tmp;
				}
				
			}else{ // bestr >= bestc
				
				
				best = bestr;
				if (r1 < r2) {
					// cumulative partial colsums
					for (j = 0; j < m; j++) {
						MAT_ELT(MS, RI[r2] ,j,n) = MAT_ELT(MT, RI[r2] ,j,n);
					}
					if (r1+1 < r2) { // not neighboring
						for (k = r2-1; k > r1; k--) {
							for (j = 0; j < m; j++) {
								MAT_ELT(MS, RI[k] ,j,n) = MAT_ELT(MT, RI[k] ,j,n) + MAT_ELT(MS, RI[k+1] ,j,n);
							}
						}
					}
					
					
					if (px == 1) {
						// changes to DC
						for (j = 0; j < m-1; j++) {
							for (j2 = j+1; j2 < m; j2++) {
								MAT_ELT(DC,j,j2,m) -= 2*( (float)MAT_ELT(MT, RI[r1] ,j2,n) * MAT_ELT(MS, RI[r1+1] ,j,n) - (float)MAT_ELT(MT, RI[r1] ,j,n) * MAT_ELT(MS, RI[r1+1] ,j2,n));
								MAT_ELT(DC,j2,j,m) = -MAT_ELT(DC,j,j2,m);
							}
						}
					}
					
					// index vector RI changes
					tmp = RI[r1];
					for (i = r1; i < r2; i++) {
						RI[i] = RI[i+1];
					}
					RI[r2] = tmp;
					
				}else {
					// cumulative partial colsums
					for (j = 0; j < m; j++) {
						MAT_ELT(MS, RI[r2] ,j,n) = MAT_ELT(MT, RI[r2] ,j,n);
					}
					if (r2+1 < r1) { // not neighboring
						for (k = r2+1; k < r1; k++) {
							for (j = 0; j < m; j++) {
								MAT_ELT(MS, RI[k] ,j,n) = MAT_ELT(MT, RI[k] ,j,n) + MAT_ELT(MS, RI[k-1] ,j,n);
							}
						}
					}
					if (px == 1) {
						// changes to DC
						for (j = 0; j < m-1; j++) {
							for (j2 = j+1; j2 < m; j2++) {
								MAT_ELT(DC,j,j2,m) -= 2*(  (float)MAT_ELT(MT, RI[r1] ,j,n) * MAT_ELT(MS, RI[r1-1] ,j2,n) - (float)MAT_ELT(MT, RI[r1] ,j2,n) * MAT_ELT(MS, RI[r1-1] ,j,n) );
								MAT_ELT(DC,j2,j,m) = -MAT_ELT(DC,j,j2,m);
							}
						}
					}
					
					// index vector RI changes
					tmp = RI[r1];
					for (i = r1; i > r2; i--) {
						RI[i] = RI[i-1];
					}
					RI[r2] = tmp;
					
					
					
				}
				
				
			}
			
			// reset
			bestc = 0;
			bestr = 0;
			
			
			crtv -= best;
			
			
			best = 0;
		}else { // best <= 0
			opt = 1;
		}
		//if( steps > 16 ){
		//	opt = 1;
		//}
		steps++;
		
	}//end opt -------------
	
	}// end for blocks
		
}// end iter -----------
	
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		REAL(out)[i] = RI[i]; //+1
	}
	for( j = 0; j < m;j++ ){
		REAL(out)[n+j] = CI[j]; //+1
	}
	
	REAL(out)[n+m] = crtv;
	free(DC2);
	free(DR2);
	free(DC);
	free(DR);
	free(MS);
	free(MT);
	//Rprintf("\nOptimization ended after %d steps with crtv=%f\n",steps,crtv);
	return out;
}


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //

// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //

SEXP quicktile(SEXP M, SEXP dims, SEXP pv){
	
	// M is byrow
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int px = INTEGER(pv)[1];
	int py = INTEGER(pv)[0];
	
	
	//Rprintf("%d\t%d\t%d\t%d\t%d\t",n,m,px,py);
	
	int i,j,i2,j2,k,tmp;
	
	
	// index vectors with initial values 0..n-1 and 0..m-1 containing the actual orders
	int RI[n];
	int CI[m];
	for (i = 0; i < n; i++) {
		RI[i] = i;
	}
	for (j = 0; j < m; j++) {
		CI[j] = j;
	}
	
	// Delta matrices for rows and columns
	// crit( i -> i2 ) - crit( i2 -> i ) if i < i2 and -(...) else
	// i.e. this is 'lost' when putting i to i2
	
	
	
	// Sum of delta matrices 
	
	
	float *DC2 = calloc(m*m,sizeof(float));
	
	float *DR2 = calloc(n*n,sizeof(float));
	float *DC = calloc(m*m,sizeof(float));
	
	float *DR = calloc(n*n,sizeof(float));
	
	int *MS = calloc(n*m,sizeof(int));	
	int *MT = calloc(n*m,sizeof(int));	
	
	
	
	if (py == 1) {
		for (i = 0; i < n; i++) {
			MAT_ELT(DR2,i,i,n) = 0;
			for (i2 = 0; i2 < n; i2++) {
				MAT_ELT(DR,i,i2,n) = 0;
			}
		}
	}
	if (px == 1) {
		for (j = 0; j < m; j++) {
			MAT_ELT(DC2,j,j,m) = 0;
			
			for (j2 = 0; j2 < m; j2++) {
				MAT_ELT(DC,j,j2,m) = 0;
			}
		}
	}
	
	//Rprintf("check02");
	// the criterion
	float crtv = 0;
	float td;
	
	//an marray pointer to M
	
	//int* mv = &INTEGER(M)[0];
	//int (*MT)[m] = (int (*)[m])mv;
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//	Rprintf("%d\t", MT[i][j]);
	//}
	//Rprintf("\n");
	//}
	
	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++) {
			MAT_ELT(MT,i,j,n) = INTEGER(M)[i+j*n];
		}
	}
	
	
	
	// --------------------- compute the initial delta values --------------------- //
	
	// cumsum by columns
	for (j = 0; j < m; j++) {
		MAT_ELT(MS,0,j,n) = MAT_ELT(MT,0,j,n);
	}
	for (i = 1; i < n; i++) {
		for (j = 0; j < m; j++) {
			MAT_ELT(MS,i,j,n) = MAT_ELT(MT,i,j,n) + MAT_ELT(MS,i-1,j,n);
		}
	}
	
	// delta C
	
	for (j = 0; j < m-1; j++) {
		for (j2 = j+1; j2 < m; j2++) {
			for (i = 1; i < n; i++) {
				td = MAT_ELT(MT,i,j,n) * MAT_ELT(MS,i-1,j2,n);
				crtv += td;
				//DC[j][j2] += ( MT[i][j]*MS[i-1][j2] - MT[i][j2]*MS[i-1][j] ); // loss
				MAT_ELT(DC,j,j2,m) += ( td - MAT_ELT(MT,i,j2,n) * MAT_ELT(MS,i-1,j,n) ); // loss
			}
			MAT_ELT(DC,j2,j,m) = - MAT_ELT(DC,j,j2,m);
		}
	}
	// also computes the initial criterion
	
	if (py == 1) {
		// cumsum by rows
		for (i = 0; i < n; i++) {
			MAT_ELT(MS,i,0,n) = MAT_ELT(MT,i,0,n);
		}
		for (i = 0; i < n; i++) {
			for (j = 1; j < m; j++) {
				MAT_ELT(MS,i,j,n) = MAT_ELT(MT,i,j,n) + MAT_ELT(MS,i,j-1,n);
			}
		}
		
		// delta R
		for (i = 0; i < n-1; i++) {
			for (i2 = i+1; i2 < n; i2++) {
				for (j = 1; j < m; j++) {
					MAT_ELT(DR,i,i2,n) += ( MAT_ELT(MT,i,j,n) * MAT_ELT(MS,i2,j-1,n) - MAT_ELT(MT,i2,j,n) * MAT_ELT(MS,i,j-1,n) ); // loss
				}
				MAT_ELT(DR,i2,i,n) = - MAT_ELT(DR,i,i2,n);
			}
		}
	}
	
	
	//Rprintf("\n\nstart crtv = %f\n\n",crtv);
	
	// --------------------- start the main algorithm --------------------- //	
	
	
	// look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to i2)
	// 'best' is the decrease in crtv, r1, r2, c1, c2 are the corresp. indices
	float bestc = 0;
	float bestr = 0;
	float best = 0;
	int r1 = 0;
	int r2 = 0;
	int c1 = 0;
	int c2 = 0;
	int opt = 0;
	int steps = 0;
	
	while( opt == 0 ){	
		
		
		//look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to k)
		
		if (py == 1) {
			for (i = 0; i < n-1; i++) {
				for (i2 = i+1; i2 < n; i2++) {
					MAT_ELT(DR2,i,i2,n) = MAT_ELT(DR2,i,i2-1,n)+ MAT_ELT(DR, RI[i] , RI[i2] ,n);
					if (MAT_ELT(DR2,i,i2,n) > bestr) {
						
						bestr = MAT_ELT(DR2,i,i2,n);
						r1 = i;
						r2 = i2;
					}
				}
				
			}
			
			for (i = 1; i < n; i++) {
				for (i2 = i-1; i2 >= 0; i2--) {
					MAT_ELT(DR2,i,i2,n) = MAT_ELT(DR2,i,i2+1,n)+  MAT_ELT(DR, RI[i2] , RI[i] ,n);
					if (MAT_ELT(DR2,i,i2,n) > bestr) {
						
						bestr = MAT_ELT(DR2,i,i2,n);
						r1 = i;
						r2 = i2;
					}
				}
				
			}
		}
		
		
		if (px == 1) {
			for (j = 0; j < m-1; j++) {
				for (j2 = j+1; j2 < m; j2++) {
					MAT_ELT(DC2,j,j2,m) = MAT_ELT(DC2,j,j2-1,m)+ MAT_ELT(DC, CI[j] , CI[j2] ,m);
					if (MAT_ELT(DC2,j,j2,m) > bestc) {
						
						bestc = MAT_ELT(DC2,j,j2,m);
						c1 = j;
						c2 = j2;
					}
				}
				
			}
			
			for (j = 1; j < m; j++) {
				for (j2 = j-1; j2 >= 0; j2--) {
					MAT_ELT(DC2,j,j2,m) = MAT_ELT(DC2,j,j2+1,m)+ MAT_ELT(DC, CI[j2] , CI[j] ,m);
					if (MAT_ELT(DC2,j,j2,m) > bestc) {
						
						bestc = MAT_ELT(DC2,j,j2,m);
						c1 = j;
						c2 = j2;
					}
				}
				
			}
		}
		
		
		
		// two columns OR two rows are exchanged. The DR/DC values will change whilst DC/DR remains unchanged.
		// the original data matrix and the delta matrices DC and DR and MS remain unchanged with repect to their row/column orders.
		// BUT the values change dependent on the permutations! Hence the sums (DR2 and DC2) use CI[j] and RI[i]
		
		//Rprintf("\n\n r1 = %d, r2 = %d, bestr = %f\n",r1,r2, bestr);
		//Rprintf("\n\n c1 = %d, c2 = %d, bestc = %f\n",c1,c2, bestc);
		
		//bestc = bestr+1;
		
		
		// reset bestc and bestr
		if( bestc > 0 || bestr > 0 ){	
			
			if (bestc > bestr) {
				best = bestc;
				
				if (c1 < c2) {
					// cumulative partial rowsums
					for (i = 0; i < n; i++) {
						MAT_ELT(MS,i, CI[c2] ,n) = MAT_ELT(MT,i, CI[c2] ,n);
					}
					if (c1+1 < c2) { // not neighboring
						for (k = c2-1; k > c1; k--) {
							for (i = 0; i < n; i++) {
								MAT_ELT(MS,i, CI[k] ,n) = MAT_ELT(MT,i, CI[k] ,n) + MAT_ELT(MS,i, CI[k+1] ,n);
							}
						}
					}
					if (py == 1) {
						// changes to DR
						for (i = 0; i < n-1; i++) {
							for (i2 = i+1; i2 < n; i2++) {
								//							Rprintf("c1 = %d, CI[c1+1] = %d\n",c1,CI[c1+1]);
								//							Rprintf("----------\n");
								//							Rprintf("MT[%d,%d] = %d, MS[%d,%d] = %d, prod = %f\n",i2,CI[c1],MT[i2][ CI[c1] ],i,CI[c1+1],MS[i][ CI[c1+1] ],2*(float)MT[i2][ CI[c1] ]*MS[i][ CI[c1+1] ]);
								//							Rprintf("----------\n");
								//							Rprintf("\nMT[%d,%d] = %d, MS[%d,%d] = %d, prod = %f\n",i,CI[c1],MT[i][ CI[c1] ],i2,CI[c1+1],MS[i2][ CI[c1+1] ],2*(float)MT[i][ CI[c1] ]*MS[i2][ CI[c1+1] ]);
								//Rprintf("----------\n");
								MAT_ELT(DR,i,i2,n) -= 2*(  (float)MAT_ELT(MT,i2, CI[c1] ,n) * MAT_ELT(MS,i, CI[c1+1] ,n) -  (float)MAT_ELT(MT,i, CI[c1] ,n) * MAT_ELT(MS,i2, CI[c1+1] ,n) );
								MAT_ELT(DR,i2,i,n) = -MAT_ELT(DR,i,i2,n);
							}
						}
					}
					
					// index vector CI changes
					tmp = CI[c1];
					for (i = c1; i < c2; i++) {
						CI[i] = CI[i+1];
					}
					CI[c2] = tmp;
					//Rprintf("------------M-------------\n");
					
					//for (i = 0; i < n; i++) {
					//	for (j = 0; j < m; j++) {
					//		Rprintf("%d\t",MT[ RI[i] ][ CI[j] ]);
					//	}
					//	Rprintf("\n");
					//}
					
					//Rprintf("----------------------------\n");
					
				}else {
					// cumulative partial rowsums
					for (i = 0; i < n; i++) {
						MAT_ELT(MS,i, CI[c2] ,n) = MAT_ELT(MT,i, CI[c2] ,n);
					}
					if (c2+1 < c1) { // not neighboring
						for (k = c2+1; k < c1; k++) {
							for (i = 0; i < n; i++) {
								MAT_ELT(MS,i, CI[k] ,n) = MAT_ELT(MT,i, CI[k] ,n) + MAT_ELT(MS,i, CI[k-1] ,n);
							}
						}
					}
					if (py == 1) {
						// changes to DR
						for (i = 0; i < n-1; i++) {
							for (i2 = i+1; i2 < n; i2++) {
								MAT_ELT(DR,i,i2,n) -= 2*(  (float)MAT_ELT(MT,i, CI[c1] ,n) * MAT_ELT(MS,i2, CI[c1-1] ,n) - (float)MAT_ELT(MT,i2, CI[c1] ,n) * MAT_ELT(MS,i, CI[c1-1] ,n) ); // neg. of c1 < c2
								MAT_ELT(DR,i2,i,n) = -MAT_ELT(DR,i,i2,n);
							}
						}
					}
					
					// index vector CI changes
					tmp = CI[c1];
					for (i = c1; i > c2; i--) {
						CI[i] = CI[i-1];
					}
					CI[c2] = tmp;
				}
				
			}else{ // bestr >= bestc
				
				
				best = bestr;
				if (r1 < r2) {
					// cumulative partial colsums
					for (j = 0; j < m; j++) {
						MAT_ELT(MS, RI[r2] ,j,n) = MAT_ELT(MT, RI[r2] ,j,n);
					}
					if (r1+1 < r2) { // not neighboring
						for (k = r2-1; k > r1; k--) {
							for (j = 0; j < m; j++) {
								MAT_ELT(MS, RI[k] ,j,n) = MAT_ELT(MT, RI[k] ,j,n) + MAT_ELT(MS, RI[k+1] ,j,n);
							}
						}
					}
					
					
					if (px == 1) {
						// changes to DC
						for (j = 0; j < m-1; j++) {
							for (j2 = j+1; j2 < m; j2++) {
								MAT_ELT(DC,j,j2,m) -= 2*( (float)MAT_ELT(MT, RI[r1] ,j2,n) * MAT_ELT(MS, RI[r1+1] ,j,n) - (float)MAT_ELT(MT, RI[r1] ,j,n) * MAT_ELT(MS, RI[r1+1] ,j2,n));
								MAT_ELT(DC,j2,j,m) = -MAT_ELT(DC,j,j2,m);
							}
						}
					}
					
					// index vector RI changes
					tmp = RI[r1];
					for (i = r1; i < r2; i++) {
						RI[i] = RI[i+1];
					}
					RI[r2] = tmp;
					
				}else {
					// cumulative partial colsums
					for (j = 0; j < m; j++) {
						MAT_ELT(MS, RI[r2] ,j,n) = MAT_ELT(MT, RI[r2] ,j,n);
					}
					if (r2+1 < r1) { // not neighboring
						for (k = r2+1; k < r1; k++) {
							for (j = 0; j < m; j++) {
								MAT_ELT(MS, RI[k] ,j,n) = MAT_ELT(MT, RI[k] ,j,n) + MAT_ELT(MS, RI[k-1] ,j,n);
							}
						}
					}
					if (px == 1) {
						// changes to DC
						for (j = 0; j < m-1; j++) {
							for (j2 = j+1; j2 < m; j2++) {
								MAT_ELT(DC,j,j2,m) -= 2*(  (float)MAT_ELT(MT, RI[r1] ,j,n) * MAT_ELT(MS, RI[r1-1] ,j2,n) - (float)MAT_ELT(MT, RI[r1] ,j2,n) * MAT_ELT(MS, RI[r1-1] ,j,n) );
								MAT_ELT(DC,j2,j,m) = -MAT_ELT(DC,j,j2,m);
							}
						}
					}
					
					// index vector RI changes
					tmp = RI[r1];
					for (i = r1; i > r2; i--) {
						RI[i] = RI[i-1];
					}
					RI[r2] = tmp;
					
					
					
				}
				
				
			}
			
			// reset
			bestc = 0;
			bestr = 0;
			
			
			crtv -= best;
			
			
			best = 0;
		}else { // best <= 0
			opt = 1;
		}
		//if( steps > 16 ){
		//	opt = 1;
		//}
		steps++;
		
	}
	
	
	
	
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		REAL(out)[i] = RI[i]; //+1
	}
	for( j = 0; j < m;j++ ){
		REAL(out)[n+j] = CI[j]; //+1
	}
	
	REAL(out)[n+m] = crtv;
	free(DC2);
	free(DR2);
	free(DC);
	free(DR);
	free(MS);
	free(MT);
	//Rprintf("\nOptimization ended after %d steps with crtv=%f\n",steps,crtv);
	return out;
}


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //


// this function will perform permutations for the categories until a local optimum is reached
SEXP quickmv(SEXP M, SEXP dims, SEXP pv , SEXP vs, SEXP minc, SEXP treevec, SEXP treelengths){
	
	// the input is as follows:
	
	// > M is the array which is to be optimized
	// > dims is a vector with the numbers of categories for each variable
	// > pv is a 0/1 vector which indicates whether or not a variable shall be reordered
	// > vs is a (currently unused) version number
	// > tree contains a (binary) tree object with 4 entries for each node:
	//  the first and the last index of the left and the right branch
	//  usually the right index of the left branch is 1 less than the left index of the right branch
	// > CM is a nd x nd matrix with optimal values for the pairwise criteria. It is used to decide which variable shall be optimized in each major step.
	// > minc (eps) is the minimal improvement which is required for a movement as a ratio of the total criterion value at the beginning of a major step.
	//  i.e. if eps = 0.01 only movements with an improvement of at least 1 percent are possible.
	
	int s,i,j,k,tmp, tmp2, d, g0, g1, g2, i1, i2, g;//r, t
	int nd = LENGTH(dims);
	
    tmp = INTEGER(treelengths)[0];
	
	int tl[nd];
    for (i = 0; i < nd; i++) {
        if (LENGTH(treelengths) > 1 ) {
            tl[i] = INTEGER(treelengths)[i];
        }else{
            tl[i] = 0;
        }
    }
        
    if( LENGTH(treelengths) > 1 ){	
		// set tree parameters and objects
        for (i = 1; i < nd; i++) {
            tmp = mymax(tmp,INTEGER(treelengths)[i]);
        }
        tmp = 3*tmp;
    }
    float *TM = calloc(tmp*nd,sizeof(int));
    //int TM[tmp][nd];
    //int *tl = &INTEGER(treelengths)[0];
	// TM is equivalent to trees, but easier to read in the code later
	// it contains the indices of the node children: min.left, min.right = max.left +1, max.right
	// TM will not change during the algorithm, the changes are made in the index vector it refers to
	
    if( LENGTH(treelengths) > 1 ){	
        k = 0;
        for (i = 0; i < tmp; i++) {
            for (j = 0; j < nd; j++) {
                MAT_ELT(TM,i,j,tmp)=0;//TM[i][j] = 0;
            }
        }
        for (d = 0; d < nd; d++) {
            for (i = 0; i < tl[d]; i++) {
                MAT_ELT(TM,3*i,d,tmp)  = INTEGER(treevec)[k + 3*i]-1;//TM[3*i][d]
                MAT_ELT(TM,3*i+1,d,tmp) = INTEGER(treevec)[k + 3*i+1]-1;//TM[3*i+1][d]
                MAT_ELT(TM,3*i+2,d,tmp) = INTEGER(treevec)[k + 3*i+2]-1;//TM[3*i+2][d]
            }
            k += 3*tl[d];
        }
    }
 	
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	//Rprintf("no. of trees = %d", nt);  
	// ml is the highest number of categories among all variables
	// biglen is the total number of categories ober all variables
	int ml = dimv[0];
	int biglen = ml;
	for(i = 1; i < nd; i++){
		biglen = biglen + dimv[i];
		if( ml < dimv[i] )
			ml = dimv[i];
	}
	// IX is a matrix with the current order of indices. Each column is for one variable
	// PX contains the current position of the original indices (only for trees)
	int IX[ml][nd];
	int PX[ml][nd];
    int cind[ml];
	// initial indices
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IX[i][k] = i;
			PX[i][k] = i;
		}
	}
	
    int *mv = &INTEGER(M)[0];
    int *IXp = &IX[0][0];
	// cp and cpr are products of the numbers of categories for the variables before and after dimension d
	// this is needed to grab the right entries from the m0 and m1 vectors
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dimv[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dimv[nd-i];
	}
	
	// the m1 array will contain the temporary (for one major step) cumulative sums over all dimensions except d
	// m0 contains the original values and is used to feed m1 in each step
	// m2 is the same as m1 except that all entries are relocated at one position further (i -> i+1)
	// all three have length mln
	
	long mln = cp[nd-1]*dimv[nd-1];
	
	float *m1 = calloc(mln,sizeof(float));
	float *m0 = calloc(mln,sizeof(float));
	float *m2 = calloc(mln,sizeof(float));
	
	
	for (i = 0; i < mln; i++) {
		m0[i] = mv[i];//INTEGER(M)[i];
	}
	for (i = 0; i < mln; i++) {
		m1[i] = mv[i];//m0[i]
		m2[i] = 0;
	}
	// TODO: compute the overall criterion alongside the computations
	// (init step and inc)
	
	float bestcrit = 0;//FLT_MAX;
	float newcrit = 0;
	float startcrit = 0;
	newcrit = mvclasscrit2(mv, IXp, dimv, nd);
	//  Rprintf("initial crit = %f",newcrit);
	int opt = 0;
	
	int  s1, s2,  j1, j2, rind, lind;//ds, s1t, s2t
	
	// DM contains the pairwise delta values
	// (i.e. the difference between the criterion when cat. 1 is before or after cat. 2)
	// DMC contains the cumulative sums over rows of DM starting from the diagonal outwards
	
	float *DM = calloc(ml*ml,sizeof(float));
	float *DMC = calloc(ml*ml,sizeof(float));
	
    // best is the best movement delta value from DMC
    // ii and jj are the corresp. indices
    // eps is the minimal value for a movement
    // eps should be either 0 or decrease with every major step
	
    float best;
    int ii, jj,  ll, lr;//gg
    float eps = 0;
	float imp;
	
	/// init done
	//Rprintf("begin optim\n");
	// d is the current dimension considered for optimization
	d = 0;
	int init = 1;
	
	// ___   MAIN STEP   ___ //
	//  |				  |  //
	//  V                 V  //   
    
	while(opt == 0){
		////
		if(INTEGER(pv)[d] > 0){
			// optimality is assumed unless an improvement can be made
			opt = 1;
			
			// current index for this dimension to work with
			for (i = 0; i < ml; i++) {
				cind[i] = i;
			}
			for (i = 0; i < ml; i++) {
				for (j = 0; j < ml; j++) {
					MAT_ELT(DM,i,j,ml) = 0;
					MAT_ELT(DMC,i,j,ml) = 0;
				}
			}
			// cumsum over all dimensions except d
			for (k = 0; k < nd; k++) {
				s1 = cp[k];
				s2 = cpr[k];
				if (d != k) {
					for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
							for (i = 1; i < dimv[k]; i++) {
								//rind = IX[i][k]*s1 + j1*s1*dimv[k] + j2;
								//lind = IX[i-1][k]*s1 + j1*s1*dimv[k] + j2;
								rind = i*s1 + j1*s1*dimv[k] + j2;
								lind = rind-s1;
								m1[ rind ] += m1[ lind ];
							}
							
						}
					}
				}
			} 
			
			//     Rprintf("test1");
			// shift the cumsums to index+1 (except dim d)
			// >> compute index increment
			int cpx[nd];
			cpx[0] = dimv[0];
			
			for (i = 1; i < nd; i++) {
				cpx[i] = cpx[i-1]*dimv[i];
			}
			int inc;
			if(d == 0){
				inc = 0;
			}else{
				inc = 1;
			}
			for (s = 0; s < nd-1; s++) {
				if( s != d-1 ){
					inc += cpx[s];
				}
			}
			//   Rprintf("inc = %d\n",inc);
			for (i = 0; i < mln-inc; i++) {
				m2[i + inc] = m1[i];
			}
			//  Rprintf("test2");
			// >> set the border (any dim expet d) to 0
			for (k = 0; k < nd; k++) {
				if(k != d){
					s1 = cp[k];
					s2 = cpr[k];
					for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
							//lind = IX[0][k]*s1 + j1*s1*dimv[k] + j2;
							lind =  j1*s1*dimv[k] + j2;
							m2[ lind ] = 0;
						}
					}
				}
			}
			
			// compute pairwise criteria from m2 and m0
			// the mv criterion is to maximize the number of 'diagonal' pairs
			// i.e. the cumsums 
			// >> first write 'i1 is before i2 to DM
			// >> (we lose that much if we put i1 after i2)
			
			s1 = cp[d];
			s2 = cpr[d];
			for (i1 = 0; i1 < dimv[d]; i1++) {
				for (i2 = 0; i2 < dimv[d]; i2++) {
					if(i1 != i2){
						for (j1 = 0; j1 < s2; j1++) {
							for (j2 = 0; j2 < s1; j2++) {
								//lind = IX[i1][d]*s1 + j1*s1*dimv[d] + j2;
								//rind = IX[i2][d]*s1 + j1*s1*dimv[d] + j2;
								lind = i1*s1 + j1*s1*dimv[d] + j2;
								rind = i2*s1 + j1*s1*dimv[d] + j2;
								MAT_ELT(DM,i1,i2,ml) += m0[ rind ] * m2[ lind ]; 
							}
						}
					}
				}
			}
			// >> now write 'put i1 after i2 difference' to DM
			
			
			for (i1 = 0; i1 < dimv[d]-1; i1++) {
				for (i2 = i1; i2 < dimv[d]; i2++) {
					if( (init == 1) & (i2 > i1) ){
						startcrit += MAT_ELT(DM,i1,i2,ml);
					}
					MAT_ELT(DM,i1,i2,ml) = MAT_ELT(DM,i2,i1,ml) - MAT_ELT(DM,i1,i2,ml);
					MAT_ELT(DM,i2,i1,ml) = -MAT_ELT(DM,i1,i2,ml);
				}
			}
			//if( (init == 1) ){
			//Rprintf("startcrit = %f\n", startcrit);
			init = 0;
			//}
			
				
			// reset m1 to m0,(m2 to 0?)
			for (i = 0; i < mln; i++) {
				m1[i] = m0[i];
				m2[i] = 0;
			}
			
			
			
			/// -----------------------------------------------------------------------		
			
			// ___   INTRA STEP 2  ___ //
			//  |				    |  //
			//  V                   V  //  
			
			// compute optimal order for veriable d from DM
			// note that DM is in the original order, DMC is in IX order
			
			// improvement must be better than eps to proceed
			//int stepcount = 0;
			imp = eps + 1;
			while( imp > eps ){
				
				if(tl[d] > 0){
					//Rprintf("imp = %f for step = %d and dim = %d\n",imp,stepcount++,d);
					imp = 0; //best = 0;
					// ---------------- TREE STEP ---------------- //
					// there is a tree for dimension d	
					
					// use IX = PX.old (IX is overwritten at the end of the tree step)
							
					for (i = 0; i < dimv[d]; i++) {
						cind[i] = IX[i][d];
						IX[i][d] = PX[i][d];
					}
					//anybetter = 0;
					
					// go through all nodes
					// both DM and the tree work on the original indices
					// i.e. we always use the same indices, but on a changed DM matrix
					// save the differences in DMC[0][ 0..tl[d]-1 ]
					for (g = 0; g < tl[d]; g++) {
						MAT_ELT(DMC,0,g,ml) = 0;
						g0 = MAT_ELT(TM,3*g+0,d,tmp);//TM[ g*3 + 0][d];
						g1 = MAT_ELT(TM,3*g+1,d,tmp);//TM[ g*3 + 1][d];
						g2 = MAT_ELT(TM,3*g+2,d,tmp);//TM[ g*3 + 2][d];
						//Rprintf("check node (%d, %d, %d):\t",PX[g0][d],PX[g1][d],PX[g2][d]);
						for (i = g0; i < g1; i++) {
							for (j = g1; j < g2+1; j++) {
								if(PX[g0][d] < PX[g1][d]){
									MAT_ELT(DMC,0,g,ml) += MAT_ELT(DM, PX[i][d] , PX[j][d] ,ml);//+= DM[ cind[i] ][ cind[j] ];
								}else{
									MAT_ELT(DMC,0,g,ml) -= MAT_ELT(DM, PX[i][d],PX[j][d] ,ml);//-= DM[ cind[i] ][ cind[j] ];
								}
							}
						}
					}
					
					// do all good steps
					// = apply change to current index and original index IX
					for (g = 0; g < tl[d]; g++) {
						// do all changes which are an improvement
						if ( MAT_ELT(DMC,0,g,ml) > 0) {
							if(MAT_ELT(DMC,0,g,ml) > best){
								imp = MAT_ELT(DMC,0,g,ml);
							}
							// ll and lr are the numbers of items in the two branches
							// both are negative iff PX[g1][d] < PX[g0][d]
							
							g0 = MAT_ELT(TM,3*g+0,d,tmp);//TM[ g*3 + 0][d];
							g1 = MAT_ELT(TM,3*g+1,d,tmp);//TM[ g*3 + 1][d];
							g2 = MAT_ELT(TM,3*g+2,d,tmp);//TM[ g*3 + 2][d];
						//	Rprintf("turning branch %d with imp = %f",g,MAT_ELT(DMC,0,g,ml));
							ll = g1-g0; //PX[g1][d]-PX[g0][d];
							lr = g2-g1+1; //PX[g2][d]-PX[g1][d]+1;
							
							bestcrit += MAT_ELT(DMC,0,g,ml);
							if (PX[g1][d] > PX[g0][d]) {
								for (s = g0; s < g1; s++) {
									PX[s][d] += lr;
								}
								for (s = g1; s < g2+1; s++) {
									PX[s][d] -= ll;
								}
							}else{
								for (s = g0; s < g1; s++) {
									PX[s][d] -= lr;
								}
								for (s = g1; s < g2+1; s++) {
									PX[s][d] += ll;
								}
							}
							
							//anybetter = 1;
							// an improvement was possible => not yet optimal
							opt = 0;
						}
					}
					//if(anybetter == 1){
					
					//rearrange m0 acc. to the new PX order
					// m1 and m0 are still the same => read from m1
					// IX is the old position, PX is the new one ( IX[i][d] == cind[ IX[i][d] ] )
					for (i = 0; i < dimv[d]; i++) {
						cind[ PX[i][d] ] = IX[i][d];
                    }
                   
					s1 = cp[d];
					s2 = cpr[d];
					
					for (i = 0; i < dimv[d]; i++) {
						for (j1 = 0; j1 < s2; j1++) {
							for (j2 = 0; j2 < s1; j2++) {
								rind = i*s1 + j1*s1*dimv[d] + j2;
								//lind = IX[i][d]*s1 + j1*s1*dimv[d] + j2;
								lind = cind[i]*s1 + j1*s1*dimv[d] + j2;
								m0[rind] = m1[lind]; 
							}
						}
					}
					for (i = 0; i < mln; i++) {
						m1[i] = m0[i];
					}
               
					// update IX (refrain from misusing IX as a temporal PX)
					for (s = 0; s < dimv[d]; s++) {
						IX[ PX[s][d] ][d] = s;
					}
					
										imp = 0;
				}else{
					// ---------------- free STEP ---------------- //
					best = 0;
					for (i = 0; i < dimv[d]-1; i++) {
						for (j = i+1; j < dimv[d]; j++) {
							MAT_ELT(DMC,i,j,ml) = MAT_ELT(DMC,i,j-1,ml) + MAT_ELT(DM,cind[i],cind[j],ml);//+ DM[ IX[i][d] ][ IX[j][d] ];
							if(MAT_ELT(DMC,i,j,ml) > best){
								ii = i;
								jj = j;
								best = MAT_ELT(DMC,i,j,ml);
							}
						}
					}
					for (i = 1; i < dimv[d]; i++) {
						for (j = i-1; j >= 0; j--) {
							MAT_ELT(DMC,i,j,ml) = MAT_ELT(DMC,i,j+1,ml) + MAT_ELT(DM,cind[j],cind[i],ml);//+ DM[ IX[j][d] ][ IX[i][d] ]; // i <-> j ?
							if(MAT_ELT(DMC,i,j,ml) > best){
								ii = i;
								jj = j;
								best = MAT_ELT(DMC,i,j,ml);
							}
						}
					}
					// do best step
					// = apply change to current index and original index IX
					if(best > eps){
						// Rprintf(" ---- ---- %d / %d is better by %f\t", ii,jj, best);
						bestcrit += best;
						tmp = IX[ii][d];
						tmp2 = cind[ii];
						if (ii < jj) {
							for (s = ii; s < jj; s++) {
								IX[s][d] = IX[s+1][d];
								cind[s] = cind[s+1];
							}
							IX[jj][d] = tmp;
							cind[jj] = tmp2;
						}
						if (ii > jj) {
							for (s = ii; s > jj; s--) {
								IX[s][d] = IX[s-1][d];
								cind[s] = cind[s-1];
							}
							IX[jj][d] = tmp;
							cind[jj] = tmp2;
						}
						// an improvement was possible => not yet optimal
						opt = 0;
					}else{
						// Rprintf("### -------- stop ------ #\n");
						
						//rearrange m0 acc. to the new IX order
						// m1 and m0 are still the same => read from m1
						s1 = cp[d];
						s2 = cpr[d];
						
						for (i = 0; i < dimv[d]; i++) {
							for (j1 = 0; j1 < s2; j1++) {
								for (j2 = 0; j2 < s1; j2++) {
									rind = i*s1 + j1*s1*dimv[d] + j2;
									//lind = IX[i][d]*s1 + j1*s1*dimv[d] + j2;
									lind = cind[i]*s1 + j1*s1*dimv[d] + j2;
									m0[rind] = m1[lind]; 
								}
							}
						}
						for (i = 0; i < mln; i++) {
							m1[i] = m0[i];
						}
					}
					imp = best;
										
				}
				
			}
			/// ''''' END INTRA STEP 2 ''''' ///
		}
		/// -----------------------------------------------------------------------				
		// decide which d is next
		d = (d+1)%nd;
		
				
	}
	/// ''''' END MAIN STEP ''''' ///
	
	
	SEXP out = allocVector(REALSXP, biglen + 1);
	
	k = 0;
	for(i = 0; i < nd; i++){
		for (j = 0; j < dimv[i]; j++) {
			REAL(out)[k] = IX[j][i];
			k++;
		}
	}
	REAL(out)[biglen] = bestcrit+newcrit;
	free(m1);
	free(m0);
	free(m2);
	free(DM);
	free(DMC);
    free(TM);
	return out;
}


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //





// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //

SEXP symmtile(SEXP M, SEXP dims, SEXP pv){
	
	// M is byrow
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	//int px = INTEGER(pv)[1];
	//int py = INTEGER(pv)[0];
	
	
	//Rprintf("%d\t%d\t%d\t%d\t%d\t",n,m,px,py);
	
	int i,j,i2,j2,k,tmp;
	
	
	// index vectors with initial values 0..n-1 and 0..m-1 containing the actual orders
	//int RI[n];
	int CI[m];
	//for (i = 0; i < n; i++) {
	//	RI[i] = i;
	//}
	for (j = 0; j < m; j++) {
		CI[j] = j;
	}
	
	// Delta matrices for rows and columns
	// crit( i -> i2 ) - crit( i2 -> i ) if i < i2 and -(...) else
	// i.e. this is 'lost' when putting i to i2
	
	
	
	// Sum of delta matrices 
	
	
	float *DC2 = calloc(m*m,sizeof(float));
	float *DC = calloc(m*m,sizeof(float));
	int *MS = calloc(n*m,sizeof(int));
	int *MT = calloc(n*m,sizeof(int));	
	
	
	
	
		//for (i = 0; i < n; i++) {
		//	MAT_ELT(DR2,i,i,n) = 0;
		//	for (i2 = 0; i2 < n; i2++) {
		//		MAT_ELT(DR,i,i2,n) = 0;
		//	}
		//}
	
	
		for (j = 0; j < m; j++) {
			MAT_ELT(DC2,j,j,m) = 0;
			
			for (j2 = 0; j2 < m; j2++) {
				MAT_ELT(DC,j,j2,m) = 0;
			}
		}
	
	
	//Rprintf("check02");
	// the criterion
	float crtv = 0;
	
	float td;
	float tmp1, tmp2;
	float crc = 0;
	
	//an marray pointer to M
	
	//int* mv = &INTEGER(M)[0];
	//int (*MT)[m] = (int (*)[m])mv;
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//	Rprintf("%d\t", MT[i][j]);
	//}
	//Rprintf("\n");
	//}
	
	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++) {
			MAT_ELT(MT,i,j,n) = INTEGER(M)[i+j*n];
		}
	}
	
	
	
	// --------------------- compute the initial delta values --------------------- //
	
	// cumsum by columns
	for (j = 0; j < m; j++) {
		MAT_ELT(MS,0,j,n) = MAT_ELT(MT,0,j,n);
	}
	for (i = 1; i < n; i++) {
		for (j = 0; j < m; j++) {
			MAT_ELT(MS,i,j,n) = MAT_ELT(MT,i,j,n) + MAT_ELT(MS,i-1,j,n);
		}
	}
	
	// delta C
	
	for (j = 0; j < m-1; j++) {
		for (j2 = j+1; j2 < m; j2++) {
			for (i = 1; i < n; i++) {
				td = MAT_ELT(MT,i,j,n) * MAT_ELT(MS,i-1,j2,n);
				crtv += td;
				//DC[j][j2] += ( MT[i][j]*MS[i-1][j2] - MT[i][j2]*MS[i-1][j] ); // loss
				MAT_ELT(DC,j,j2,m) += ( td - MAT_ELT(MT,i,j2,n) * MAT_ELT(MS,i-1,j,n) ); // loss
			}
			MAT_ELT(DC,j,j2,m) = 2*(MAT_ELT(DC,j,j2,m) + MAT_ELT(MT,j,j,n) * MAT_ELT(MT,j2,j2,n) - MAT_ELT(MT,j,j2,n) * MAT_ELT(MT,j2,j,n));
           // MAT_ELT(DC,j,j2,m) = 2*( MAT_ELT(DC,j,j2,m));
			MAT_ELT(DC,j2,j,m) = - MAT_ELT(DC,j,j2,m);
		}
	}
   // Rprintf("\n");
	//for (i=0; i<m; i++) {
      //  for (j=0; j<m; j++) {
        //    Rprintf("%f\t",MAT_ELT(DC,i,j,m));
        //}
        //Rprintf("\n");
    //}
    //Rprintf("\n");
	
	
//	Rprintf("\n\nstart crtv = %f\n\n",crtv);
//Rprintf("Symmetric");
	// --------------------- start the main algorithm --------------------- //	
	
	
	// look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to i2)
	// 'best' is the decrease in crtv, r1, r2, c1, c2 are the corresp. indices
	float bestc = 0;
	//float bestr = 0;
	float best = 0;
	int c1 = 0;
	int c2 = 0;
	int opt = 0;
	int steps = 0;
	
	while( opt == 0 ){	
		
		
			//if (px == 1) {
			for (j = 0; j < m-1; j++) {
				for (j2 = j+1; j2 < m; j2++) {
					tmp1 = 0;
					tmp2 = 0;
					if (j2-j > 1) {
						for (i = j+1; i < j2; i++) {
							tmp1 += MAT_ELT(MT, CI[i] , CI[j] ,m);
							tmp2 += MAT_ELT(MT, CI[i] , CI[j2] ,m);
						}
						tmp1 = tmp1 *  MAT_ELT(MT, CI[j] , CI[j2] ,m);
						tmp2 = tmp2 * MAT_ELT(MT, CI[j] , CI[j] ,m);
					}
					MAT_ELT(DC2,j,j2,m) = MAT_ELT(DC2,j,j2-1,m)+ MAT_ELT(DC, CI[j] , CI[j2] ,m) - 4*tmp1 + 4*tmp2;
					if (MAT_ELT(DC2,j,j2,m) > bestc) {
						bestc = MAT_ELT(DC2,j,j2,m);
						c1 = j;
						c2 = j2;
					}
				}
				
			}
			
			for (j = 1; j < m; j++) {
				for (j2 = j-1; j2 >= 0; j2--) {
					tmp1 = 0;
					tmp2 = 0;
					if (j-j2 > 1) {
						for (i = j2+1; i < j; i++) {
							tmp1 += MAT_ELT(MT, CI[i] , CI[j2] ,m);
							tmp2 += MAT_ELT(MT, CI[i] , CI[j] ,m);
						}
						tmp1 = tmp1 * MAT_ELT(MT, CI[j] , CI[j] ,m);
						tmp2 = tmp2 * MAT_ELT(MT, CI[j2] , CI[j] ,m);
					}
					
					
					MAT_ELT(DC2,j,j2,m) = MAT_ELT(DC2,j,j2+1,m)+ MAT_ELT(DC, CI[j2] , CI[j] ,m) + 4*tmp1 - 4*tmp2;
					if (MAT_ELT(DC2,j,j2,m) > bestc) {
						bestc = MAT_ELT(DC2,j,j2,m);
						c1 = j;
						c2 = j2;
					}
				}
				
			}
		//}
		
		

		// reset bestc and bestr
		//if( bestc > 0 || bestr > 0 ){	
			if( bestc > 0 ){
			
				best = bestc;
				
				if (c1 < c2) {
					// cumulative partial rowsums
					for (i = 0; i < n; i++) {
						MAT_ELT(MS,i, CI[c2] ,n) = MAT_ELT(MT,i, CI[c2] ,n);
					}
					if (c1+1 < c2) { // not neighboring
						for (k = c2-1; k > c1; k--) {
							for (i = 0; i < n; i++) {
								MAT_ELT(MS,i, CI[k] ,n) = MAT_ELT(MT,i, CI[k] ,n) + MAT_ELT(MS,i, CI[k+1] ,n);
							}
						}
					}
					
					// symmetric matrix: 
										
					// changes to DC = DR+DC
					for (i = 0; i < n-1; i++) {
						for (i2 = i+1; i2 < n; i2++) {
                            
                            MAT_ELT(DC,i,i2,n) -= 4*(  (float)MAT_ELT(MT,i2, CI[c1] ,n) * MAT_ELT(MS,i, CI[c1+1] ,n) -  (float)MAT_ELT(MT,i, CI[c1] ,n) * MAT_ELT(MS,i2, CI[c1+1] ,n) );
                            if(i == CI[c1]){
                                for (j= c1+1; j <= c2; j++) {
                                    if(CI[j] == i2){
                                        MAT_ELT(DC,i,i2,n) += 4*( (float)MAT_ELT(MT,i2, i ,n) * MAT_ELT(MT, i, i2 ,n)- (float)MAT_ELT(MT,i, i ,n) * MAT_ELT(MT, i2, i2 ,n));
                                    }
                                }
                            }
                           
                            MAT_ELT(DC,i2,i,n) = -MAT_ELT(DC,i,i2,n);
						}
					}
                    
					
									
					// index vector CI changes
					tmp = CI[c1];
					for (i = c1; i < c2; i++) {
						CI[i] = CI[i+1];
					}
					CI[c2] = tmp;
					//Rprintf("------------M-------------\n");
					
					//for (i = 0; i < n; i++) {
					//	for (j = 0; j < m; j++) {
					//		Rprintf("%d\t",MT[ RI[i] ][ CI[j] ]);
					//	}
					//	Rprintf("\n");
					//}
					
					//Rprintf("----------------------------\n");
					
				}
				if(c1 > c2) {
					// cumulative partial rowsums
					for (i = 0; i < n; i++) {
						MAT_ELT(MS,i, CI[c2] ,n) = MAT_ELT(MT,i, CI[c2] ,n);
					}
					if (c2+1 < c1) { // not neighboring
						for (k = c2+1; k < c1; k++) {
							for (i = 0; i < n; i++) {
								MAT_ELT(MS,i, CI[k] ,n) = MAT_ELT(MT,i, CI[k] ,n) + MAT_ELT(MS,i, CI[k-1] ,n);
							}
						}
					}
						
					for (i = 0; i < n-1; i++) {
						for (i2 = i+1; i2 < n; i2++) {
							MAT_ELT(DC,i,i2,n) -= 4*(  (float)MAT_ELT(MT,i, CI[c1] ,n) * MAT_ELT(MS,i2, CI[c1-1] ,n) - (float)MAT_ELT(MT,i2, CI[c1] ,n) * MAT_ELT(MS,i, CI[c1-1] ,n) ); // neg. of c1 < c2
							if(i2 == CI[c1]){
                                for (j= c1-1; j >= c2; j--) {
                                    if(CI[j] == i){
                                        MAT_ELT(DC,i,i2,n) += 4*( (float)MAT_ELT(MT,i2, i ,n) * MAT_ELT(MT, i, i2 ,n)- (float)MAT_ELT(MT,i, i ,n) * MAT_ELT(MT, i2, i2 ,n));
                                    }
                                }
                            }
                           
                            MAT_ELT(DC,i2,i,n) = -MAT_ELT(DC,i,i2,n);
						}
					}
                   
					
					// index vector CI changes
					tmp = CI[c1];
					for (i = c1; i > c2; i--) {
						CI[i] = CI[i-1];
					}
					CI[c2] = tmp;
				}
				
			
			
			// reset
			bestc = 0;
			
			
			//Rprintf("\n#c1 = %d,c2 = %d\n",c1,c2);
			crtv -= best;
				crc -= best;
  
			
			best = 0;
		}else { // best <= 0
			opt = 1;
		}
if( steps > 2*m ){
	opt = 1;
}
		steps++;
        

        
			
	}
	
	
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		REAL(out)[i] = CI[i]; //+1
	}
	for( j = 0; j < m;j++ ){
		REAL(out)[n+j] = CI[j]; //+1
	}
	
	REAL(out)[n+m] = crtv;
	free(DC2);
	
	free(DC);
	
	free(MS);
	free(MT);
	//Rprintf("\nOptimization ended after %d steps with crtv=%f\n",steps,crtv);
	return out;
}


// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------------------------- //


