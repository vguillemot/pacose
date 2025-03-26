#include <R.h>
#include <Rmath.h>

void wermuthC(int *p, int *ninteract, double *delta, double *error, int *iter, int *maxiter,
                    int ind[*ninteract][2], double icovx[*p][*p], double result[*p][*p])
{
 double sii,sjj,sij,d,errortemp,coef ;
 int i,j,k,l,u,i1,i2 ;
 double temp[*p][*p] ;
 for (i1 = 0; i1 < *p; i1++) for (i2 = 0; i2 < *p; i2++)  temp[i1][i2] = icovx[i1][i2];
 while( (*iter < *maxiter) && (*error > *delta) ){
   *iter = *iter+1 ;
   for (u = 0; u < *ninteract; u++) {
        i = fmin2(ind[u][0],ind[u][1]);
        j = fmax2(ind[u][0],ind[u][1]);
        sii = temp[i][i] ;
        sjj = temp[j][j] ;
        sij = temp[i][j] ;
        d   = sii*sjj-sij*sij ;
        result[i][j] = result[j][i] = 0 ;
        result[i][i] = d/sjj ;
        result[j][j] = d/sii ;
        for (k=0; k < (*p) ; k++) {
          if (k != i && k != j) {
            result[k][i] = result[i][k] = temp[i][k] - sij*temp[j][k]/sjj;
            result[k][j] = result[j][k] = temp[j][k] - sij*temp[i][k]/sii;
            for (l=0; l < (*p) ; l++) {
              if (l != i && l != j && l <= k) result[k][l] = result[l][k] = temp[k][l]
                                                  - (sij/d)*( (temp[i][k])*(temp[j][l] - sij*temp[i][l]/sjj) +
                                                              (temp[j][k])*(temp[i][l] - sij*temp[j][l]/sjj) ) ;
            }
          }
        }
   for (i1 = 0; i1 < *p; i1++) for (i2 = 0; i2 < *p; i2++) temp[i1][i2] = result[i1][i2];
   }
   errortemp = 0;
   for (u = 0; u < *ninteract; u++) {
     i = ind[u][0] ;
     j = ind[u][1] ;
     coef =  fmax2(-result[i][j],result[i][j]) ;
     //printf ("Coef: %d , %d , %f \n", i,j,coef);
     if (coef > errortemp) errortemp = coef;
   }
   *error = errortemp ;
 }
}

