#if !defined(MNEWT_H__INCLUDED_)
#define MNEWT_H__INCLUDED_

int mnewtc(int ntrial, double x[], int n, double tolf, double minx[], double maxx[],
           void (*usrfun)(double *, double **, double *, double *), double usrfunpar[], double *berr, int nlos, double tollos, RANDOMDEF *rd);

#endif
