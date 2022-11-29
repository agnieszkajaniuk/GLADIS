#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "konfig.h"
#include "dysk_zwm.h"
#include "steady.h"
#include "nrutil.h"
#include "mnewt.h"
#include "random.h"
#include "boundary.h"

//static int  clearusrfun;

/* void usrfun_mnewt_scurve(double x[], double **lalpha, double *lbetha, double *usrfunpar); */
void moja(double xpom[], double *usrfunpar, int mrozno);



/**************************************************************************/
double scurve (FORONERAD *stab)
/*calculates stability curve for a given radius in Schwarzschild units*/
/**************************************************************************/
{
    int i;
    double Tc, rhoc;
    double Omega, /*Omegamag,*/ fluxtot, etaj;
    double mdottab[MROZ], Mdtab[MROZ];
    double /*a5, a6,*/ Hd, Te;
    double  xpom[4];
    //    double ro/*, T*/;
    //double tolerancje[5];
    //FILE *plik1;
    double prom;
    //, Prommag;
    double xi, Alpha2;

    prom = stab->r;

    //printf("Calculating scurve at prom = %E\n", prom);

    //    Prommag = stab->Rmag;


    //    clearusrfun = 0;                           //Obliczanie od nowa danych dla usrfun

    Tc = 3.5e4*pow((1.e8/konf.Mass*Msun),0.25);      /*initial values of disc temperature and density*/
    rhoc = (5.e-8)*(1.e8/konf.Mass*Msun);

    xpom[1] = rhoc;
    xpom[2] = Tc;

    //    Alpha2=Alpha;
    if(konf.viscosity) {
        xi=4./3.*Sigmabol/Crad*pow(xpom[2], 3.)*Mprot/Kbol/xpom[1];
        //       Alpha2=Flicker(prom, 0.0)*(1+xi/konf.xi0)/(1+pow((xi/konf.xi0), 2));
        Alpha2 = konf.Alpha*(1+xi/konf.xi0)/(1+pow((xi/konf.xi0), 2));
    } else {

        //        Alpha2=Flicker(prom, 0.0);
        Alpha2 = konf.Alpha;
    }

    Omega = pow(Graw*konf.Mass/pow(prom*Rschw,3.0), 0.5);


    //a5 = Kbol/Mprot/Omega/Omega;
    //a6 = 4.0*Sigmabol/Crad/3.0/Omega/Omega;

    stab->kep = Omega;

    stab->taudyn = konf.Kappa_dyn/Omega;

    //    stab->taudynmag = Kappa_dyn/Omegamag;

    mdottab[0] = -3.0;
    mdottab[MROZ-1] = 3.0;

    for (i=1; i<=MROZ-2; i++) {
        mdottab[i]=mdottab[i-1]+(mdottab[MROZ-1]-mdottab[0])/MROZ;

    }


    for (i=0; i < MROZ; i++)
        //for (i=MROZ-1; i>=0; i--)

    {

        Mdtab[i]=pow(10, mdottab[i])*16.0*Ledd/Crad/Crad;

        fluxtot = 3.0/8.0/Pi*Graw*konf.Mass*Mdtab[i]/pow(prom*Rschw, 3.0)*(1.0-sqrt(3.0/prom));

        //(1-pow((3./prom),1.5)*(prom-1.)/2.);

        fluxtot=fluxtot/konf.C4;                         /*to agree with vert. str. coefficients*/



        etaj=1.-1./(1.+konf.Dzet*pow((pow(10,mdottab[i])/2.),2));



        //tolerancje[1] = xpom[1]/1.0E+6;   //Tolerancja dla ro
        //tolerancje[2] = xpom[2]/1.0E+6;   //Tolerancja dla T

        double usrfunpar[6];

        usrfunpar[1] = fluxtot;
        usrfunpar[2] = prom;
        usrfunpar[3] = etaj;
        usrfunpar[4] = Alpha2;


        //printf("mdot %E\n", mdottab[i]);

        moja(xpom, usrfunpar, i);

        //printf("Po mnewt %E %E berr=%E\n", xpom[1], xpom[2], berr);
        //getc(stdin);

        //2008.01.30 Hd = sqrt((a5*xpom[2]+a6*pow(xpom[2], 4.0)/(xpom[1]))/konf.C3);
        Hd = xpom[3];

        //x=mdottab[i];
        //y=log10(xpom[1]*Hd);

        /*modified viscosity law*/
        if(prom == konf.Rmax) {
            //Alpha2=Flicker(prom, 0.0);
            Alpha2 = konf.Alpha;
        } else {
            if(konf.viscosity) {
                xi=4./3.*Sigmabol/Crad*pow(xpom[2], 3.)*Mprot/Kbol/xpom[1];
                //                Alpha2=Flicker(prom, 0.0)*(1+xi/konf.xi0)/(1+pow((xi/konf.xi0), 2));
                Alpha2 = konf.Alpha*(1+xi/konf.xi0)/(1+pow((xi/konf.xi0), 2));
            } else {

                //                Alpha2=Flicker(prom, 0.0);
                Alpha2 = konf.Alpha;
            }

            if(konf.hotvis) {
                Alpha2=konf.Alpha0*pow((Hd/prom/Rschw), -0.5);

                //                if(prom < 10.0)
                //                {printf("prom=%lf, T=%lf, alpha=%lf\n", prom, xpom[2], Alpha2);}
            }
            else {
                //                Alpha2=Flicker(prom, 0.0);
                Alpha2 = konf.Alpha;
            }



        }
        stab->sigma[i] = xpom[1]*Hd;
        stab->temper[i] = xpom[2];
        stab->h[i] = Hd;
        //          stab->xvr[i] = 0.5*Mdtab[i]/Pi/prom;

        stab->etajet[i] = etaj;
        stab->Mdotstab[i] = pow(10, mdottab[i]);

        if(konf.heating) {
            Te=pow(fluxtot/Sigmabol*(1-konf.Xicor)*(1.0-etaj), 0.25);
        } else {
            Te=pow(fluxtot/Sigmabol, 0.25);
        }

        stab->Tef[i] = Te;

        //if(prom > 4.0 && prom < 5.0)
        //  { printf("prom = %E, Te=%E, etaj=%E\n", prom, Te, etaj);
        //getc(stdin);
        //  }
        //          fprintf(plik1,"%10f  %10f  %10f \n",curtm.Teff, xpom[2]/pow(0.75*(xpom[1]*Hd*Kappa+0.66), 0.25), xpom[2]);

    }
    //       fclose(plik1);

    return 0;

}
/************************************************************************/
static void usrfun_mnewt_scurve_ipsqrt(double x[], double **lalpha, double *lbetha, double *pars)
/************************************************************************/
{
    double Hd, beta, P, Pgas, ro, T;
    double F1, F1_eps;
    double eps = x[1]*1.0e-6;

    double etaj    = pars[0];
    double fluxtot = pars[1];
    double Omega   = pars[2];
    double a7      = pars[3];
    double a4      = pars[4];
    double Alpha2  = pars[5];



    Hd = x[1];
    beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;

    // heating propto sqrt(Pgas*Ptot): sqrt(beta)

    P = fluxtot/konf.C1/Hd/Alpha2/a7/sqrt(beta);
    Pgas = beta*P;
    ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
    T = Pgas*Mprot/Kbol/ro;
    F1 = (P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0));

    ///////
    Hd = x[1] + eps;
    beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;

    // heating propto sqrt(Pgas*Ptot): sqrt(beta)

    P = fluxtot/konf.C1/Hd/Alpha2/a7/sqrt(beta);
    Pgas = beta*P;
    ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
    T = Pgas*Mprot/Kbol/ro;
    F1_eps = (P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0));

    lalpha[1][1] = (F1_eps - F1)/eps;

    lbetha[1] = -F1;

}

/************************************************************************/
static void usrfun_mnewt_scurve_iptot(double x[], double **lalpha, double *lbetha, double *pars)
/************************************************************************/
{
    double Hd, beta, P, Pgas, ro, T;
    double F1, F1_eps;
    double eps = x[1]*1.0e-6;

    double etaj    = pars[0];
    double fluxtot = pars[1];
    double Omega   = pars[2];
    double a7      = pars[3];
    double a4      = pars[4];
    double Alpha2  = pars[5];
    double a8      = pars[6];
    double a3      = pars[7];



    Hd = x[1];
    P = fluxtot/konf.C1/Hd/Alpha2/a7;
    ro = P*a8/konf.C3/Hd/Hd;
    T = fluxtot*((1.0-konf.Xicor)*(1.0-etaj) - Hd*Hd*a4)/konf.C2/a3;
    T = T *ro*Hd;
    T = pow(T, 0.25);
    F1 = P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0);


    ///////
    Hd = x[1] + eps;
    P = fluxtot/konf.C1/Hd/Alpha2/a7;
    ro = P*a8/konf.C3/Hd/Hd;
    T = fluxtot*((1.0-konf.Xicor)*(1.0-etaj) - Hd*Hd*a4)/konf.C2/a3;
    T = T *ro*Hd;
    T = pow(T, 0.25);
    F1_eps = P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0);

    lalpha[1][1] = (F1_eps - F1)/eps;

    lbetha[1] = -F1;

}

/************************************************************************/
static void usrfun_mnewt_scurve_ipgas(double x[], double **lalpha, double *lbetha, double *pars)
/************************************************************************/
{
    double Hd, beta, P, Pgas, ro, T;
    double F1, F1_eps;
    double eps = x[1]*1.0e-6;

    double etaj    = pars[0];
    double fluxtot = pars[1];
    double Omega   = pars[2];
    double a7      = pars[3];
    double a4      = pars[4];
    double Alpha2  = pars[5];
    double a8      = pars[6];
    double a3      = pars[7];



    Hd = x[1];
    Pgas = fluxtot/konf.C1/Hd/Alpha2/a7;
    beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;
    ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
    T = Pgas*Mprot/Kbol/ro;
    P = Pgas/beta;
    F1 = (P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0))/P;


    ///////
    Hd = x[1] + eps;
    Pgas = fluxtot/konf.C1/Hd/Alpha2/a7;
    beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;
    ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
    T = Pgas*Mprot/Kbol/ro;
    P = Pgas/beta;
    F1_eps = (P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0))/P;

    lalpha[1][1] = (F1_eps - F1)/eps;

    lbetha[1] = -F1;

}

/************************************************************************/
void moja(double xpom[], double *usrfunpar, int mrozno)
/*Vertically averaged equations to calculate the stability curve*/
/************************************************************************/
{
    double /*a1, a2,*/ a3, a4, /*a5, a6,*/ a7, a8;
    double Omega;

    double Hd, F1, F1min, Hdbest;
    double ro;
    double T;
    double P, Pgas, beta;
    //    double a1_new, a2_new;
    double Hmin, Hmax, dHl;
    int i, imax;

    double fluxtot = usrfunpar[1];
    double prom = usrfunpar[2];
    double etaj = usrfunpar[3];
    double Alpha2 = usrfunpar[4];

    const int Wesja_Bozenki = 0;

    imax = 20000;

    ro = xpom[1];
    T = xpom[2];
    Hdbest = 1e35;


    //    printf("x[1] = %G, x[2] = %G\n", x[1], x[2]);

    //    if (clearusrfun == 0) {

    //        clearusrfun = 1;

    Omega = pow(Graw*konf.Mass/pow(prom*Rschw,3.0), 0.5);

    //a1 = 3.0/2.0*Omega*Kbol/Mprot;
    //a2 = 4.0/2.0*Omega*Sigmabol/Crad;

    a3 = 4.0*Sigmabol/3.0/konf.Kappa;
    a4 = 0.5*4.0/3.0*konf.Qadv/(1.0-sqrt(3.0/prom))/prom/prom/Rschw/Rschw;
    //(1-pow((3./prom),1.5)*(prom-1.)/2.)/prom/prom/Rschw/Rschw;
    //a5 = Kbol/Mprot/Omega/Omega;
    //a6 = 4.0*Sigmabol/Crad/3.0/Omega/Omega;
    a7 = 3.0/2.0*Omega;
    a8 = 1.0/Omega/Omega;
    //    }


    //    a1_new=a1*Alpha2;
    //a2_new=a2*Alpha2;

    if(konf.ipsqrt) {

        //	     Hmin = fluxtot*(1.0-konf.Xicor)*(1.0-etaj)*konf.Kappa/konf.C2/konf.C3/Crad/Omega/Omega;

        //Hmin corrected: 01/29/07

        double A = (1.0 - konf.Xicor)*(1.0 - etaj);
        double C = fluxtot/konf.C2/konf.C3/Crad/Omega/Omega*konf.Kappa;
        Hmin = ( -1.0 + sqrt(1.0 + 4.0*A*C*C*a4))/(2*a4*C);
        Hmax = pow((1.0-konf.Xicor)*(1.0-etaj)/a4, 0.5);


        if(Wesja_Bozenki)
        {
            dHl = (log10(Hmax)-log10(Hmin))/(imax - 1.0);

            F1min = 1.e35;

            for(i=1; i<=imax; i++) {
                Hd = pow(10,log10(Hmin)+dHl*(i-1));
                beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;

                // heating propto sqrt(Pgas*Ptot): sqrt(beta)

                P = fluxtot/konf.C1/Hd/Alpha2/a7/sqrt(beta);
                Pgas = beta*P;
                ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
                T = Pgas*Mprot/Kbol/ro;
                F1 = (P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0))/P;
                if(F1*F1<F1min*F1min) {
                    F1min = F1;
                    Hdbest = Hd;
                }
            }
            Hd = Hdbest;
        }
        else
        {
            double pars[6];
            double berr = 1.0e100;

            pars[0] = etaj;
            pars[1] = fluxtot;
            pars[2] = Omega;
            pars[3] = a7;
            pars[4] = a4;
            pars[5] = Alpha2;

            double xint[2];
            double minx[2], maxx[2];

            dHl = (log10(Hmax)-log10(Hmin))/(imax - 1.0);

            xint[1] = pow(10,log10(Hmin)+dHl*(imax/2-1));

            minx[1] = Hmin;
            maxx[1] = Hmax;

            RANDOMDEF rd;
            setRandomTable(&rd, 12345*mrozno, 2345*mrozno, 34533*mrozno, 84533*mrozno);

            mnewtc(200, xint, 1, 1.0e-7, minx, maxx, usrfun_mnewt_scurve_ipsqrt, pars, &berr, 200, 1.0e-6, &rd);

            Hd = xint[1];
        }

        beta = 1.0 - ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;

        P = fluxtot/konf.C1/Hd/Alpha2/a7/sqrt(beta);
        Pgas = beta*P;
        ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
        T = Pgas*Mprot/Kbol/ro;
    } else {

        if(konf.iptot) {

            Hmin = (prom - 3.0)*Rschw/10000.;
            Hmax = pow((1.0-konf.Xicor)*(1.0-etaj)/a4, 0.5);

            if(Wesja_Bozenki)
            {
                dHl = (log10(Hmax)-log10(Hmin))/(imax - 1.0);
                F1min = 1.e35;

                for(i=1; i<=imax; i++) {
                    Hd = pow(10,log10(Hmin)+dHl*(i-1));
                    P = fluxtot/konf.C1/Hd/Alpha2/a7;
                    ro = P*a8/konf.C3/Hd/Hd;
                    T = fluxtot*((1.0-konf.Xicor)*(1.0-etaj) - Hd*Hd*a4)/konf.C2/a3;
                    T = T *ro*Hd;
                    T = pow(T, 0.25);
                    F1 = P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0);
                    if(F1*F1<F1min*F1min) {
                        F1min = F1;
                        Hdbest = Hd;
                    }
                }
                Hd = Hdbest;
            }
            else
            {

                double pars[8];
                double berr = 1.0e100;

                pars[0] = etaj;
                pars[1] = fluxtot;
                pars[2] = Omega;
                pars[3] = a7;
                pars[4] = a4;
                pars[5] = Alpha2;
                pars[6] = a8;
                pars[7] = a3;

                double xint[2];
                double minx[2], maxx[2];

                dHl = (log10(Hmax)-log10(Hmin))/(imax - 1.0);

                xint[1] = pow(10,log10(Hmin)+dHl*(imax/2-1));

                minx[1] = Hmin;
                maxx[1] = Hmax;

                RANDOMDEF rd;
                setRandomTable(&rd, 12345*mrozno, 2345*mrozno, 34533*mrozno, 84533*mrozno);

                mnewtc(200, xint, 1, 1.0e-7, minx, maxx, usrfun_mnewt_scurve_iptot, pars, &berr, 200, 1.0e-6, &rd);

                Hd = xint[1];

            }


            P = fluxtot/konf.C1/Hd/Alpha2/a7;
            ro = P*a8/konf.C3/Hd/Hd;
            T = fluxtot*((1.0-konf.Xicor)*(1.0-etaj) - Hd*Hd*a4)/konf.C2/a3;
            T = T *ro*Hd;
            T = pow(T, 0.25);
        } else {

            //       Hmin = fluxtot*(1.0-konf.Xicor)*(1.0-etaj)*konf.Kappa/konf.C2/konf.C3/Crad/Omega/Omega;

            //Hmin corrected: 01/29/07
            double A = (1.0 - konf.Xicor)*(1.0 - etaj);
            double C = fluxtot/konf.C2/konf.C3/Crad/Omega/Omega*konf.Kappa;
            Hmin = ( -1.0 + sqrt(1.0 + 4.0*A*C*C*a4))/(2*a4*C);
            Hmax = pow((1.0-konf.Xicor)*(1.0-etaj)/a4, 0.5);

            if(Wesja_Bozenki)
            {

                dHl = (log10(Hmax)-log10(Hmin))/(imax - 1.0);
                F1min = 1.e35;

                for(i=1; i<=imax; i++) {
                    Hd = pow(10,log10(Hmin)+dHl*(i-1));
                    Pgas = fluxtot/konf.C1/Hd/Alpha2/a7;
                    beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;
                    ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
                    T = Pgas*Mprot/Kbol/ro;
                    P = Pgas/beta;
                    F1 = (P - Kbol/Mprot*ro*T - 4.0*Sigmabol/Crad/3.0*pow(T, 4.0))/P;
                    if(F1*F1<F1min*F1min) {
                        F1min = F1;
                        Hdbest = Hd;
                    }
                    //	 printf("%4d  Pg=%e  b=%e  ro=%e  T=%e  F1=%e   F1min=%e \n",i,Pgas,beta,ro,T,F1,F1min);
                }
                Hd = Hdbest;
            }
            else
            {
                double pars[8];
                double berr = 1.0e100;

                pars[0] = etaj;
                pars[1] = fluxtot;
                pars[2] = Omega;
                pars[3] = a7;
                pars[4] = a4;
                pars[5] = Alpha2;
                pars[6] = a8;
                pars[7] = a3;

                double xint[2];
                double minx[2], maxx[2];

                dHl = (log10(Hmax)-log10(Hmin))/(imax - 1.0);

                xint[1] = pow(10,log10(Hmin)+dHl*(imax/2-1));

                minx[1] = Hmin;
                maxx[1] = Hmax;

                RANDOMDEF rd;
                setRandomTable(&rd, 12345*mrozno, 2345*mrozno, 34533*mrozno, 84533*mrozno);

                mnewtc(200, xint, 1, 1.0e-7, minx, maxx, usrfun_mnewt_scurve_ipgas, pars, &berr, 200, 1.0e-6, &rd);

                Hd = xint[1];
            }


            Pgas = fluxtot/konf.C1/Hd/Alpha2/a7;
            beta = 1.0- ((1.0-konf.Xicor)*(1.0-etaj)-Hd*Hd*a4)*fluxtot/konf.C2/konf.C3/Crad/Omega/Omega/Hd*konf.Kappa;
            ro = Pgas/beta/konf.C3/Hd/Hd/Omega/Omega;
            T = Pgas*Mprot/Kbol/ro;
            // printf("r=%e  Pg=%e  b=%.2lf  ro=%e  T=%e   Hd=%e  flux=%e \n",prom,Pgas,beta,ro,T,Hd,fluxtot);

        }
    }
    xpom[1] = ro;
    xpom[2] = T;
    xpom[3] = Hd;
    //    exit(1);
    //    if(prom>3.5) {exit(1);}
    return ;
}
