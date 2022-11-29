#ifdef MPI_USED
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


#include "konfig.h"
#include "dysk_zwm.h"
#include "boundary.h"
#include "mdothold.h"

extern double Markoff(double urand, RANDOMDEF *rd);


static int PeriodicBoundaryMdot(double *mpar, SOLINTIME *curtm, double time_);
static int TabularizedBoundaryMdot(double *mpar, SOLINTIME *curtm, double time_);
static void ReCalcBoundary(SOLINTIME *curtm, double mdote, double time_);

static double toIndicesInTab(double md);



//Tabularized mdot
static MdotHold *llm0dat = NULL;

/**************************************************************/
void InitOuterBoundary(void)
/**************************************************************/
{

    if(konf.external_mdot == 2)
        llm0dat = new MdotHold("llm0.dat");

}

/**************************************************************/
void DeInitOuterBoundary(void)
/**************************************************************/
{
    if(llm0dat)
    {
        delete llm0dat;
        llm0dat = NULL;
    }
}


/**************************************************************/
int OuterBoundary(double *mpar, SOLINTIME *curtm, double time_)
/**************************************************************/
{
    switch(konf.external_mdot) {

    case 0: //constant boundary contidions, do nothing
        return 0;

    case 1: //periodic boundary contidions for mdot
        return PeriodicBoundaryMdot(mpar, curtm, time_);

    case 2: //tabularized boundary contidions for mdot
        return TabularizedBoundaryMdot(mpar, curtm, time_);
    }

    return 0;
}


/**************************************************************/
static int PeriodicBoundaryMdot(double *mpar, SOLINTIME *curtm, double time_)
/**************************************************************/
{
    double mdote;

    //	mdote = konf.Mdot*(1 - sin(2.0*Pi/konf.Tper*time_)/2.0);

    mdote = sin(2.0*Pi/konf.Tper*time_);
    mdote = pow(10, mdote);
    mdote = konf.Mdot*mdote;

    //    *mpar = mdote*Ledd/Crad/Crad*16.0;
    *mpar = mdote*Msun/Year;

    ReCalcBoundary(curtm, mdote, time_);

    return 0;
}



/**************************************************************/
static int TabularizedBoundaryMdot(double *mpar, SOLINTIME *curtm, double time_)
/**************************************************************/
{
    double mdote;

    mdote = llm0dat->findMdot(time_);       //search for mdot in llm0.dat

    //	 printf("mdote= %E\n", mdote);

    *mpar = mdote*Ledd/Crad/Crad*16.0;

    ReCalcBoundary(curtm, mdote, time_);

    return 0;
}
/**************************************************************/
static void ReCalcBoundary(SOLINTIME *curtm, double mdote, double time_)
/**************************************************************/
{

    double ji;
    ji = toIndicesInTab(mdote);

    //	printf("ji= %E\n", ji);

    curtm->sig_ster[PROZ-1]= splineInterp(gr.tab[PROZ-1].Mdotstab, gr.tab[PROZ-1].sigma, MROZ, ji);
    curtm->temp_ster[PROZ-1]= splineInterp(gr.tab[PROZ-1].Mdotstab, gr.tab[PROZ-1].temper, MROZ, ji);

    curtm->s[PROZ-1]=curtm->sig_ster[PROZ-1]*gr.y[PROZ-1];

    curtm->hd_ster[PROZ-1]=  splineInterp(gr.tab[PROZ-1].Mdotstab, gr.tab[PROZ-1].h, MROZ, ji);;
    curtm->nu_ster[PROZ-1]=0.66*Flicker(PROZ-1,curtm->urandalfa[PROZ-1], gr.tab[PROZ-1].r*Rschw, time_, curtm->hd_ster[PROZ-1])/gr.tab[PROZ-1].kep;
    curtm->nu_ster[PROZ-1]=curtm->nu_ster[PROZ-1]*(curtm->hd_ster[PROZ-1]/curtm->sig_ster[PROZ-1]*4/3*Sigmabol/Crad
                           *pow(curtm->temp_ster[PROZ-1], 4.0) +
                           Kbol/Mprot*curtm->temp_ster[PROZ-1]);
    curtm->snu[PROZ-1]=curtm->s[PROZ-1]*curtm->nu_ster[PROZ-1];

    curtm->Ptot[PROZ-1]=Kbol/Mprot*curtm->sig_ster[PROZ-1]/curtm->hd_ster[PROZ-1]
                        *curtm->temp_ster[PROZ-1]+4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[PROZ-1],4.0);

    double ksi=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[PROZ-1],3.0)*Mprot/Kbol*curtm->hd_ster[PROZ-1]/curtm->sig_ster[PROZ-1];
    curtm->beta[PROZ-1]=1.0/(1.0+ksi);
    curtm->hd_old[PROZ-1]=curtm->hd_ster[PROZ-1];

    curtm->Pgas[PROZ-1]=Kbol/Mprot*curtm->sig_ster[PROZ-1]/curtm->hd_ster[PROZ-1]
                        *curtm->temp_ster[PROZ-1];
    curtm->Prad[PROZ-1]=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[PROZ-1],4.0);
    curtm->vr_dysk[PROZ-1]= 6.0/curtm->s[PROZ-1]/gr.y[PROZ-1]*(curtm->snu[PROZ-1]-curtm->snu[PROZ-2])/gr.dy[PROZ-1];

    if(konf.wymiana_masy) {

        curtm->sig_cor[PROZ-1]=1.0/konf.Kappa;
        //curtm->sig_cor[PROZ-1]=0.1/konf.Kappa; // glebokosc opt. korony na r_out;
        curtm->h_cor[PROZ-1] = gr.tab[PROZ-1].r*Rschw;
        curtm->t_cor[PROZ-1] = Graw*Mprot/Kbol*konf.Mass/gr.tab[PROZ-1].r/Rschw;
        curtm->s_cor[PROZ-1] = curtm->sig_cor[PROZ-1]*gr.y[PROZ-1];
        curtm->nu_cor[PROZ-1]= 0.66*konf.AlphaCor/gr.tab[PROZ-1].kep*Kbol/Mprot*curtm->t_cor[PROZ-1];
        curtm->snu_cor[PROZ-1] = curtm->s_cor[PROZ-1]*curtm->nu_cor[PROZ-1];
        curtm->lasttimedynprom[PROZ-1] = 0.0;
    }


    curtm->flux[PROZ-1] = 0.0;
    curtm->flux2[PROZ-1] = 0.0;
    curtm->flux3[PROZ-1] = 0.0;
    curtm->fluxcor[PROZ-1]=0.0;

    curtm->urand[PROZ-1] = 0.0;
    curtm->Alpha22[PROZ-1]=Flicker(PROZ-1,curtm->urandalfa[PROZ-1], gr.tab[PROZ-1].r*Rschw, time_,curtm->hd_ster[PROZ-1]);

}


/*****************************************************************/
static double toIndicesInTab(double md)
/*****************************************************************/
{

    double ret;

    ret = md;

    if(ret < 1.0e-3)
    {
        printf("toIndicesInTab: mdot lower than limit: %E\n", ret);
        ret = 1.0e-3;
    }

    if(ret > 1.0e+3)
    {
        printf("toIndicesInTab: mdot grater than limit: %E\n", ret);
        ret = 1.0e+3;
    }

    return ret;
}

/******************************************************************************/
double CalcUrandAlfa(int jr, double urandalfa, double time_, RANDOMDEF *rd, SOLINTIME *newtm )
/******************************************************************************/
{
    int jm;
    double usum;

    usum = 0.0;

    
    int lastjm = gr.tabmagalfaRange[jr].last;
    
    int ile = lastjm - gr.tabmagalfaRange[jr].first + 1;
    
    if(konf.flicker_visc) {
    	
      double P2 = 1.0/newtm->Alpha22[jr]*pow(gr.tab[jr].r*Rschw/newtm->hd_ster[jr], 2.0);
	    
		for (jm=gr.tabmagalfaRange[jr].first; jm <= lastjm; ++jm) {

			gr.tabmagalfa[jm].taudynmag = gr.tabmagalfa[jm].P1*P2;

			if (time_ >= gr.tabmagalfa[jm].lasttimedyn + gr.tabmagalfa[jm].taudynmag) {

				gr.tabmagalfa[jm].lasttimedyn = time_;

				gr.tabmagalfa[jm].urand = Markoff(gr.tabmagalfa[jm].urand, rd);
			}
			
			usum += gr.tabmagalfa[jm].urand;
		}
    	
    }
	else {
		for (jm=gr.tabmagalfaRange[jr].first; jm <= lastjm; ++jm) {
			
			if (time_ >= gr.tabmagalfa[jm].lasttimedyn + gr.tabmagalfa[jm].taudynmag) {

				gr.tabmagalfa[jm].lasttimedyn = time_;

				gr.tabmagalfa[jm].urand = Markoff(gr.tabmagalfa[jm].urand, rd);
			}
			
			usum += gr.tabmagalfa[jm].urand;
		}
	}


    if(ile>0)
    {
        usum = usum/ile;
        
        if(usum != urandalfa)
        {
	  //    	    printf("CalcUrandAlpha changed jr=%2d time=%E newval=%E, oldval=%E\n", jr, time_, usum, urandalfa);
        }
        
        return usum;
    }

    return urandalfa;
}




/****************************************************************************/
double Flicker(int jr, double urandalfa, double prom, double time_, double H)
/****************************************************************************/
{

    //  double bpar;
    double Beta;
    double AlphaFlick;
    //  double tauFlick;

    if(konf.alfa_flickering) {

        //  printf("prom*Rschw=%E, urand=%lf\n", prom, urandalfa);

        // bpar = 0.01;
        // bpar = 0.1;
        /*bpar = 0.99;*/

        // tauFlick = 1.0/gr.tab[jr].kep*pow((prom/hd_ster), 2.0);

        //  AlphaFlick = konf.Alpha*(1.0 + Apar/tauFlick*exp(-1.0*time_/tauFlick)*pow(prom, bpar));

        //alpha(r,t) = alph0 (1+ A/tau exp(-t/tau) r^(b))
        //tau = 1/(alpha0 Omega_K(r))(r/H)^2


      /* Correction to mimic the changing H-grid while disk in outburst*/

      double H0 = gr.tab[jr].Rmagalfa/0.98; //as read from the file rsig.dat
      //      double H = curtm->hd_ster[jr];

      Beta = konf.bpar*urandalfa/sqrt(H0/H);
      //      printf("t=%E, jr=%d, H0=%E, H=%E, Beta=%E\n", time_, jr, H0, H, Beta);

        AlphaFlick = konf.Alpha*(1.0 + Beta);
        
        const double AlphaMin = 1.0e-10;

        if (AlphaFlick < AlphaMin)
        	AlphaFlick = AlphaMin;        
        //  AlphaFlick = konf.Alpha;

        //    printf("t=%E, jr=%d, urandalfa=%E, alpha=%E\n", time_, jr, urandalfa, AlphaFlick);
        return AlphaFlick;

    }
    else {

        return konf.Alpha;
    }

}

/****************************************************************************/
void Boundary(SOLINTIME *curtm)
/****************************************************************************/
{
    if(rank == 0)
    {
        const int jr = 0;
        curtm->temp_ster[jr] = curtm->temp_ster[1];

//#if 0
        curtm->sig_ster[jr]  = curtm->s[1]/gr.y[1];
        curtm->s[jr]         = curtm->sig_ster[jr]*gr.y[jr];
        curtm->Ptot[jr]      = curtm->Ptot[1];

        curtm->hd_ster[jr]=curtm->Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/curtm->sig_ster[jr]/sqrt(konf.C3);

        double dtemp = Kbol/Mprot*curtm->sig_ster[jr]/curtm->hd_ster[jr]*curtm->temp_ster[jr];

        curtm->Ptot[jr]= dtemp + 4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr], 4.0);
        curtm->Pgas[jr]= dtemp;
        curtm->Prad[jr]=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr], 4.0);
        curtm->hd_ster[jr]=curtm->Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/curtm->sig_ster[jr]/sqrt(konf.C3);

        if(konf.ipsqrt) {
            curtm->nu_ster[jr]=sqrt(curtm->Pgas[jr]*curtm->Ptot[jr])*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*curtm->hd_ster[jr]/curtm->sig_ster[jr];

        } else {
            if (konf.iptot) {
                curtm->nu_ster[jr]=curtm->Ptot[jr]*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*curtm->hd_ster[jr]/curtm->sig_ster[jr];
            } else {
                curtm->nu_ster[jr]=curtm->Pgas[jr]*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*curtm->hd_ster[jr]/curtm->sig_ster[jr];
            }

        }

        curtm->snu[jr]=curtm->s[jr]*curtm->nu_ster[jr];
//#endif

        if(konf.wymiana_masy) 
        {
            curtm->sig_cor[jr] = curtm->s_cor[1]/gr.y[1];
            curtm->s_cor[jr]   = curtm->sig_cor[jr]*gr.y[jr];
            curtm->nu_cor[jr]  = 0.66*konf.AlphaCor/gr.tab[jr].kep;
            curtm->nu_cor[jr] *= Kbol/Mprot*curtm->t_cor[jr];
            curtm->snu_cor[jr] = curtm->s_cor[jr]*curtm->nu_cor[jr];
        }
    }


    if(rank == size - 1)
    {
        const int jr = PROZ-1;
        curtm->temp_ster[jr] = curtm->temp_ster[PROZ-2];

//#if 0
        curtm->sig_ster[jr]  = curtm->s[PROZ-2]/gr.y[PROZ-2];
        curtm->s[jr]         = curtm->sig_ster[jr]*gr.y[jr];

        curtm->Ptot[jr]      = curtm->Ptot[PROZ-2];

        curtm->hd_ster[jr]=curtm->Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/curtm->sig_ster[jr]/sqrt(konf.C3);

        double dtemp = Kbol/Mprot*curtm->sig_ster[jr]/curtm->hd_ster[jr]*curtm->temp_ster[jr];

        curtm->Ptot[jr]= dtemp + 4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr], 4.0);
        curtm->Pgas[jr]= dtemp;
        curtm->Prad[jr]=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr], 4.0);
        curtm->hd_ster[jr]=curtm->Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/curtm->sig_ster[jr]/sqrt(konf.C3);

        if(konf.ipsqrt) {
            curtm->nu_ster[jr]=sqrt(curtm->Pgas[jr]*curtm->Ptot[jr])*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*curtm->hd_ster[jr]/curtm->sig_ster[jr];

        } else {
            if (konf.iptot) {
                curtm->nu_ster[jr]=curtm->Ptot[jr]*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*curtm->hd_ster[jr]/curtm->sig_ster[jr];
            } else {
                curtm->nu_ster[jr]=curtm->Pgas[jr]*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*curtm->hd_ster[jr]/curtm->sig_ster[jr];
            }

        }

        curtm->snu[jr]=curtm->s[jr]*curtm->nu_ster[jr];
//#endif

        if(konf.wymiana_masy)
        {
            curtm->sig_cor[jr] = curtm->s_cor[PROZ-2]/gr.y[PROZ-2];
            curtm->s_cor[jr]   = curtm->sig_cor[jr]*gr.y[jr];
            curtm->nu_cor[jr]  = 0.66*konf.AlphaCor/gr.tab[jr].kep;
            curtm->nu_cor[jr] *= Kbol/Mprot*curtm->t_cor[jr];
            curtm->snu_cor[jr] = curtm->s_cor[jr]*curtm->nu_cor[jr];
        }
    }
}
