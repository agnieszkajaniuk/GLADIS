#include <math.h>
#include <stdio.h>

#include "nrutil.h"
#include "random.h"

#define FREERETURN {free_matrix(alpha,1,n,1,n);free_vector(bet,1,n);\
free_ivector(indx,1,n); return 1;}

double **matrix(int,int,int,int);
int ludcmp(double **a,int n,int *indx,double *d);
void lubksb(double **a,int n,int *indx, double *b);





int mnewtc(int ntrial, double x[], int n, double tolf, double minx[], double maxx[],
           void (*usrfun)(double *, double **, double *, double *), double usrfunpar[], double *berr, int nlos, double tollos, RANDOMDEF *rd) {
    int k, i, *indx, l;
    double errf = 1.0E100;
    double d,*bet;
    double **alpha;
    double bestval[n+1];
    double besterr = 1.0E100;


    int retval = 0;


    indx=ivector(1,n);
    bet=vector(1,n);
    alpha=matrix(1,n,1,n);

    for(l = 0; l < nlos; ++l) {
        for (k=1; k<=ntrial; k++) {

            usrfun(x,alpha, bet, usrfunpar);
            errf=0.0;
            for (i=1; i<=n; i++)
                errf += fabs(bet[i]);

            //printf("mnet x1=%E x2=%E errf=%E\n", x[1], x[2], errf);
            //getc(stdin);


            if (errf <= tolf) {
                *berr = errf;
                retval = 1;
                goto koniec;
            }

            if(errf < besterr) {
                if(k > 1)
                    retval = 2;

                besterr = errf;
                for (i=1; i<=n; i++)
                    bestval[i] = x[i];
            }

            if(ludcmp(alpha,n,indx,&d))
                lubksb(alpha,n,indx,bet);

            for (i=1; i<=n; i++) {
                for(int dv = 1; dv <= 4; dv *= 2) {

                    double newval = x[i] + bet[i]/dv;

                    if(newval >= minx[i] && newval <= maxx[i]) {
                        x[i] = newval;
                        break;
                    }
                }
            }
        }

        if(besterr <= tollos)
            break;

        for (i=1; i<=n; i++) {
//            double newval = bestval[i]*pow(10.0, (0.5-(double)rand()/RAND_MAX)*1.5);
//            if(newval >= minx[i] && newval <= maxx[i])
//                x[i] = newval;

            x[i] = UniformRandomRangeLog(rd, minx[i], maxx[i]);

        }

    }


    for (i=1; i<=n; i++)
        x[i] = bestval[i];

    *berr = besterr;


koniec:
    free_matrix(alpha,1,n,1,n);
    free_vector(bet,1,n);
    free_ivector(indx,1,n);
    return retval;
}



#undef FREERETURN







