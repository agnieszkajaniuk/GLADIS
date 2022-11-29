#ifdef MPI_USED
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "konfig.h"
#include "dysk_zwm.h"
#include "nsolve.h"

/* maximum fractional increase in timestep per timestep */
#define SAFE	(1.3)
#define dtMin   (1.e-5)

int EstimateStep(double *dt, SOLINTIME *curtm, SOLINTIME *newtm);


SOLVESTRUCT solve;

/**************************************************************/
int Euler(double *dt, SOLINTIME *curtm, SOLINTIME *newtm)
/*
Dziala ma currtm, wypluwa newtm
*/
/**************************************************************/
{
    int ret = 1;

    DERIVAT ktemp;

    derivatives(*dt, curtm, &ktemp);
    calculations(*dt, curtm, newtm, &ktemp, *dt, 1);

    if(!konf.krok_staly)
        ret = EstimateStep(dt, curtm, newtm);

    return ret;
}

/**************************************************************/
int Heun(double *dt, SOLINTIME *curtm, SOLINTIME *newtm)
/*
Dziala ma currtm, wypluwa newtm
*/
/**************************************************************/
{
    int ret = 1, jr;

    DERIVAT k1d;
    DERIVAT k2d;
    DERIVAT ksumd;

    derivatives(*dt, curtm, &k1d);
    calculations(*dt, curtm, newtm, &k1d, *dt, 0);

    derivatives(*dt/2.0, newtm, &k2d);

    for (jr = konf.ji; jr <= konf.jo; jr++) {
        ksumd.dsdt[jr] = k1d.dsdt[jr] + k2d.dsdt[jr];
        ksumd.dtempdt[jr] = k1d.dtempdt[jr] + k2d.dtempdt[jr];

        if(konf.wymiana_masy) {
            ksumd.dsdt_cor[jr] = k1d.dsdt_cor[jr] + k2d.dsdt_cor[jr];
        }

    }

    calculations(*dt, curtm, newtm, &ksumd, *dt/2.0, 1);


    if(!konf.krok_staly)
        ret = EstimateStep(dt, curtm, newtm);

    return ret;
}



/**************************************************************/
int  RungeKutta(double *dt, SOLINTIME *curtm, SOLINTIME *newtm)
//Fourth order Runge-Kutta method, takes curtm, gives newtm
/**************************************************************/
{
    register int jr;
    int          ret = 1;

    static DERIVAT k1d;
    static DERIVAT k2d;
    static DERIVAT k3d;
    static DERIVAT k4d;
    static DERIVAT ksumd;

    derivatives(*dt, curtm, &k1d); //Obliczamy k1 (pochodna w punkcie)
    calculations(*dt/2.0, curtm, newtm, &k1d, *dt/2.0, 0);


    derivatives(*dt, newtm, &k2d);
    calculations(*dt/2.0, curtm, newtm, &k2d, *dt/2.0, 0);


    derivatives(*dt, newtm, &k3d);
    calculations(*dt, curtm, newtm, &k3d, *dt, 0);


    derivatives(*dt, newtm, &k4d);


    //Obliczamy sume wspolczynnikow
    for (jr = konf.ji; jr <= konf.jo; jr++) {
        ksumd.dsdt[jr] = (k1d.dsdt[jr] + k2d.dsdt[jr]*2.0
                          + k3d.dsdt[jr]*2.0 + k4d.dsdt[jr])/6.0;
        ksumd.dtempdt[jr] = (k1d.dtempdt[jr] + k2d.dtempdt[jr]*2.0
                             + k3d.dtempdt[jr]*2.0 + k4d.dtempdt[jr])/6.0;

        if(konf.wymiana_masy) {
            ksumd.dsdt_cor[jr] = (k1d.dsdt_cor[jr] + k2d.dsdt_cor[jr]*2.0
                                  + k3d.dsdt_cor[jr]*2.0 + k4d.dsdt_cor[jr])/6.0;
        }

    }


    calculations(*dt, curtm, newtm, &ksumd, *dt, 1);

    if(!konf.krok_staly)
        ret = EstimateStep(dt, curtm, newtm);

    return ret;
}




/**************************************************************/
int PredictorCorrector(double *dt, int ntime, SOLINTIME *curtm, SOLINTIME *newtm)
/**************************************************************/
{

    SOLINTIME pktm;

    DERIVAT dersum;
    DERIVAT derakt;
    DERIVAT dertemp;

    int jr;
    int    ret;
    double tdt;
    double newdt;
    double accuracy;


    if (ntime < 10000) {  //Pierwsze N punktow liczymy rk4


        tdt = *dt;

        derivatives(*dt, newtm, &dertemp);

        ret = RungeKutta(&tdt, curtm, newtm);

        //        ret = Euler(&tdt, curtm, newtm);

        if(!ret)
            return ret;


        memcpy(&solve.derm3, &solve.derm2, sizeof(DERIVAT));
        memcpy(&solve.derm2, &solve.derm1, sizeof(DERIVAT));
        memcpy(&solve.derm1, &dertemp, sizeof(DERIVAT));

        memcpy(&solve.currm2, &solve.currm1, sizeof(SOLINTIME));
        memcpy(&solve.currm1, curtm, sizeof(SOLINTIME));

        return 1;
    }


    derivatives(*dt, curtm, &derakt);

    for (jr = konf.ji; jr <= konf.jo; jr++) {
        dersum.dsdt[jr] = (55.0*derakt.dsdt[jr] - 59.0*solve.derm1.dsdt[jr]
                           + 37.0*solve.derm2.dsdt[jr]- 9.0*solve.derm3.dsdt[jr])/24.0;

        dersum.dtempdt[jr] = (55.0*derakt.dtempdt[jr] - 59.0*solve.derm1.dtempdt[jr]
                              + 37.0*solve.derm2.dtempdt[jr] - 9.0*solve.derm3.dtempdt[jr])/24.0;

        if(konf.wymiana_masy) {
            dersum.dsdt_cor[jr] = (55.0*derakt.dsdt_cor[jr] - 59.0*solve.derm1.dsdt_cor[jr]
                                   + 37.0*solve.derm2.dsdt_cor[jr] - 9.0*solve.derm3.dsdt_cor[jr])/24.0;
        }
    }

    calculations(*dt, curtm, &pktm, &dersum, *dt, 0);
    derivatives(*dt, &pktm, &dertemp);

    for (jr = konf.ji; jr <= konf.jo; jr++) {
        dersum.dsdt[jr] =  (9.0*dertemp.dsdt[jr] + 19.0*derakt.dsdt[jr]
                            - 5.0*solve.derm1.dsdt[jr] + solve.derm2.dsdt[jr])/24.0;
        dersum.dtempdt[jr] = (9.0*dertemp.dtempdt[jr] + 19.0*derakt.dtempdt[jr]
                              - 5.0*solve.derm1.dtempdt[jr] + solve.derm2.dtempdt[jr])/24.0;

        if(konf.wymiana_masy) {
            dersum.dsdt_cor[jr] =  (9.0*dertemp.dsdt_cor[jr] + 19.0*derakt.dsdt_cor[jr]
                                    - 5.0*solve.derm1.dsdt_cor[jr] + solve.derm2.dsdt_cor[jr])/24.0;
        }
    }

    calculations(*dt, curtm, newtm, &dersum, *dt, 1);


    if(konf.krok_staly)
        goto wynikok;

    newdt = *dt;
    ret = EstimateStep(&newdt, &pktm, newtm);


    if(ret == 0)
        newdt = *dt/2.0;


    if(*dt > newdt)  //Deviding the time-step
    {

        *dt = *dt/2.0;

        //        printf("Dividing the time-step dt = %E %d\n", *dt, ntime);

        memcpy(&solve.derm2, &solve.derm1, sizeof(DERIVAT));

        //        for (jr = 0; jr < PROZ; jr++)
        //            temptm.hd_old[jr] = temptm.hd_ster[jr];  //Przyblizenie liniowe
        //  curtm.hd_old[jr] = (curtm.hd_ster[jr]-curtm.hd_old[jr])/2.0;

        tdt = *dt;
        RungeKutta(&tdt, &solve.currm2, newtm);
        derivatives(*dt, newtm, &solve.derm3);

        memcpy(&solve.currm2, &solve.currm1, sizeof(SOLINTIME));

        //        for (jr = 0; jr < PROZ; jr++)
        //            solve.currm2.hd_old[jr]=newtm->hd_ster[jr];

        tdt = *dt;
        RungeKutta(&tdt, &solve.currm2, &solve.currm1);
        derivatives(*dt, &solve.currm1, &solve.derm1);

        return 0;  //ponownie liczymy krok
    }

    else if(*dt < newdt) //Multipying the time-step
    {
        SOLINTIME temptm;

        //memcpy(&curtm, &starttm, sizeof(SOLINTIME));
        tdt = (*dt)*2.0;

        ret = RungeKutta(&tdt, newtm, &temptm);
        //        ret = Euler(&tdt, newtm, &temptm);

        if(ret == 0) {
            if(rank == 0)
                printf("Cannot increase the time-step dt = %E\n", (*dt)*2.0);
            goto wynikok;

        }

        memcpy(&solve.derm2, &solve.derm1, sizeof(DERIVAT));

        derivatives(*dt, newtm, &solve.derm1);

        memcpy(&solve.currm2, &solve.currm1, sizeof(SOLINTIME));

        memcpy(&solve.currm1, newtm, sizeof(SOLINTIME));

        memcpy(newtm, &temptm, sizeof(SOLINTIME));

        ptime = ptime + *dt;

        *dt = (*dt)*2.0;

        //        printf("Multiplying the time-step dt = %E %d\n", *dt, ntime);


        //        for (jr = 0; jr < PROZ; jr++)
        //        {
        //            currm0.hd_old[jr]=solve.currm1.hd_ster[jr];
        //        }
    }
    else {

wynikok:

        memcpy(&solve.derm3, &solve.derm2, sizeof(DERIVAT));
        memcpy(&solve.derm2, &solve.derm1, sizeof(DERIVAT));
        memcpy(&solve.derm1, &derakt, sizeof(DERIVAT));

        memcpy(&solve.currm2, &solve.currm1, sizeof(SOLINTIME));
        memcpy(&solve.currm1, curtm, sizeof(SOLINTIME));
    }


    return 1;

}


int EstimateStep(double *dt, SOLINTIME *curtm, SOLINTIME *newtm) {

    double ndt = 1.e99;
    const double courant = 0.8;

    for(int j = konf.ji; j <= konf.jo; ++j) {
        double cval = courant*gr.dy[j]/fabs(newtm->vr_dysk[j]);

        if (cval < ndt)
          ndt = cval;

		if(konf.wymiana_masy && (!konf.PBK || (konf.PBK && (ptime > konf.PBK_time)))) {
            cval = courant*gr.dy[j]/fabs(newtm->vr_cor[j]);

            if (cval < ndt)
                ndt = cval;
        }
    }

#ifdef MPI_USED
    //Pobieramy najmniejsze ndt wszystkich watkow
    if(size > 1) {
        double  min_ndt;
        MPI_Allreduce( &ndt, &min_ndt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

        ndt = min_ndt;
    }

#endif
        switch(konf.nsolve_method) {

        case 0:
        case 1:
        case 3:
            if(ndt > (*dt)*SAFE)
              *dt = (*dt)*SAFE;
            else if(ndt < dtMin)
              *dt = dtMin;
            else
              *dt = ndt;

            break;
        case 2:
            if(ndt >= (*dt)*2.*SAFE)
                *dt *= 2.0;
            else if(ndt <= *dt)
                *dt /= 2.0;
            break;
        }

    return 1;
}

/**************************************************************/
void DerivDr1(double *f, double *vr, double *dfdr)
/**************************************************************/
{
    int j;
   // const int typ = 3;
      const int typ = 0;

    static int CaclCoef1 = 0;
    static double coef1[PROZ];
    static double coef2[PROZ];
    static double coef3[PROZ];
    static double coef4[PROZ];


    if(!CaclCoef1) {

        for(j = konf.ji; j <= konf.jo; ++j)
        {
            coef1[j] = (gr.y[j] - gr.y[j+1])/(gr.y[j-1] - gr.y[j])/(gr.y[j-1] - gr.y[j+1]);
            coef2[j] = (2.0*gr.y[j] - gr.y[j-1] - gr.y[j+1])/(gr.y[j] - gr.y[j-1])/(gr.y[j] - gr.y[j+1]);
            coef3[j] = (gr.y[j] - gr.y[j-1])/(gr.y[j+1] - gr.y[j-1])/(gr.y[j+1] - gr.y[j]);
        }
        for(j = konf.ji; j <= konf.jo+1; ++j)
        {
            coef4[j] = 1.0/gr.dy[j];
        }

        CaclCoef1 = 1;
    }


    switch(typ) {
    case 0: //central
        for(j = konf.ji; j <= konf.jo; ++j)
            dfdr[j] = (f[j+1] - f[j-1])/2.0*coef4[j];
        break;
    case 1: //First Derivative for unequally spaced data, Three-Point Formula
        for(j = konf.ji; j <= konf.jo; ++j)
            dfdr[j] = f[j-1]*coef1[j] + f[j]*coef2[j] + f[j+1]*coef3[j];
        break;
    case 3:   //asymetrycznie, zeleznie od vr

        for(j = konf.ji; j <= konf.jo; ++j) {

            if(vr[j] >= 0.0)
                dfdr[j] = (f[j] - f[j-1])*coef4[j];
            else
                dfdr[j] = (f[j+1] - f[j])*coef4[j+1];
        }
        break;

    }

}


/**************************************************************/
void DerivDr2(double *f, double *dfdr)
/**************************************************************/
{
    int j;
    const int typ = 0;

    static int CaclCoef2 = 0;
    static double coef21[PROZ];
    static double coef22[PROZ];
    static double coef23[PROZ];
    static double coef2[PROZ];


    if(!CaclCoef2) {

        for(j = konf.ji; j <= konf.jo; ++j) {

            coef21[j] = 2.0/(gr.y[j-1] - gr.y[j])/(gr.y[j-1] - gr.y[j+1]);
            coef22[j] = 2.0/(gr.y[j] - gr.y[j-1])/(gr.y[j] - gr.y[j+1]);
            coef23[j] = 2.0/(gr.y[j+1] - gr.y[j-1])/(gr.y[j+1] - gr.y[j]);
//printf("COEF j=%d 21=%E 22=%E 23=%E\n", j, coef21[j], coef22[j], coef23[j] );
			coef2[j] = 1.0/gr.dy[j]/gr.dy[j];
        }

        CaclCoef2 = 1;
    }


    switch(typ) {
    case 0: //central
        for(j = konf.ji; j <= konf.jo; ++j)
            dfdr[j] = (f[j+1] + f[j-1] - 2.0*f[j])*coef2[j];
        break;
    case 1: //Second Derivative for unequally spaced data
        for(j = konf.ji; j <= konf.jo; ++j)
        {
            dfdr[j] = f[j-1]*coef21[j] + f[j]*coef22[j] + f[j+1]*coef23[j];
//		      printf("DFDR j=%d, %E %E %E %E %E %E %E\n",dfdr[j],f[j-1],coef21[j], f[j], coef22[j], f[j+1], coef23[j]);
        }

        break;
    }

}
