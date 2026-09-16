#include <R.h>
#include <math.h>
#include <string.h>

void wermuthC(int *p, int *ninteract, double *delta, double *error, int *iter, int *maxiter,
                      int ind[*ninteract][2], double icovx[*p][*p], double result[*p][*p])
{
  double sii,sjj,sij,d,inv_sii,inv_sjj,sij_over_d,errortemp,coef ;
  int i,j,k,l,u ;
  const int pp = *p ;
  double temp[pp][pp] ;
  memcpy(temp, icovx, (size_t)pp * (size_t)pp * sizeof(double));
 while( (*iter < *maxiter) && (*error > *delta) ){
   *iter = *iter+1 ;
   for (u = 0; u < *ninteract; u++) {
         i = ind[u][0] < ind[u][1] ? ind[u][0] : ind[u][1];
         j = ind[u][0] < ind[u][1] ? ind[u][1] : ind[u][0];
         double *temp_i = temp[i] ;
         double *temp_j = temp[j] ;
         double *result_i = result[i] ;
         double *result_j = result[j] ;
         sii = temp_i[i] ;
         sjj = temp_j[j] ;
         sij = temp_i[j] ;
        d   = sii*sjj-sij*sij ;
         inv_sii = 1.0 / sii ;
         inv_sjj = 1.0 / sjj ;
         sij_over_d = sij / d ;
        result_i[j] = result_j[i] = 0 ;
         result_i[i] = d * inv_sjj ;
         result_j[j] = d * inv_sii ;
        for (k=0; k < pp ; k++) {
          if (k != i && k != j) {
            /* hoisted out of the l-loop: constant for the whole row k */
            double tik = temp_i[k] ;
            double tjk = temp_j[k] ;
            double *temp_k = temp[k] ;
            double *result_k = result[k] ;
            result_k[i] = result_i[k] = tik - sij*tjk*inv_sjj;
            result_k[j] = result_j[k] = tjk - sij*tik*inv_sii;
            for (l=0; l <= k ; l++) {
              if (l != i && l != j) {
                double value = temp_k[l] - sij_over_d * ( tik*(temp_j[l] - sij*temp_i[l]*inv_sjj) +
                                                           tjk*(temp_i[l] - sij*temp_j[l]*inv_sjj) ) ;
                result_k[l] = value ;
                result[l][k] = value ;
              }
            }
          }
        }
  memcpy(temp, result, (size_t)pp * (size_t)pp * sizeof(double));
   }
   errortemp = 0;
   for (u = 0; u < *ninteract; u++) {
    i = ind[u][0] ;
    j = ind[u][1] ;
    coef = fabs(result[i][j]) ;
     if (coef > errortemp) errortemp = coef;
   }
  *error = errortemp ;
 }
}

