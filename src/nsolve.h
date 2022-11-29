#if !defined(NSOLVE_H__INCLUDED_)
#define NSOLVE_H__INCLUDED_


typedef struct {

    SOLINTIME currm1;
    SOLINTIME currm2;

    DERIVAT derm1;
    DERIVAT derm2;
    DERIVAT derm3;

}
SOLVESTRUCT;

extern SOLVESTRUCT solve;

int Euler(double *dt, SOLINTIME *curtm, SOLINTIME *newtm);
int RungeKutta(double *dt, SOLINTIME *curtm, SOLINTIME *newtm);
int PredictorCorrector(double *dt, int ntime, SOLINTIME *curtm, SOLINTIME *newtm);
int Heun(double *dt, SOLINTIME *curtm, SOLINTIME *newtm);

void DerivDr1(double *f, double *vr, double *dfdr);
void DerivDr2(double *f, double *dfdr);


#endif
