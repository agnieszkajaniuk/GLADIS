#if !defined(DYSK_ZWM_H__INCLUDED_)
#define DYSK_ZWM_H__INCLUDED_

#include "random.h"

typedef struct {
    double sigma[MROZ];
    double temper[MROZ];
    double h[MROZ];
    double Tef[MROZ];
    double r;
    double Rmag;
    double Rmagalfa;
    double kep;
    double taudyn;
    //double taudynmag;
    double etajet[MROZ];
    double Mdotstab[MROZ];
    double rr;

}
FORONERAD;


typedef struct {

    double sig_ster[PROZ], temp_ster[PROZ], hd_ster[PROZ], hd_old[PROZ];
    double s[PROZ], snu[PROZ];
    double nu_ster[PROZ], signu[PROZ];
    double Teff[PROZ], flux[PROZ], logflux, lum;
    double Ptot[PROZ], beta[PROZ];
    double Pgas[PROZ], Prad[PROZ], Qplus[PROZ], Qminus[PROZ];
    double vr_dysk[PROZ];  //Predkosc radialna w dysku


    // WYMIANAMASY

    /*Parmetry zwiazane z wymiana masy z korona*/
    double sig_cor[PROZ];   //Gestosc powierzchniowa w koronie
    double s_cor[PROZ];     //Pomocnicza
    double snu_cor[PROZ];   //S*nu
    double nu_cor[PROZ];    //Lepkosc korony

    double h_cor[PROZ];  //Grubosc korony rowna promieniowi
    double mz[PROZ];     //Tempo parowania
    double Bz[PROZ];     //Pole magnetyczne w kier. z
    double Valvf[PROZ];  //Predkosc Alvfena w dysku
    double BzMax[PROZ];

    double vr_cor[PROZ];  //Predkosc radialna w koronie
    double t_cor[PROZ];  //Temperatura korony rowna temperaturze wirialnej
    //    double fluxcor[PROZ];
    //    double lumcor;


    // !WYMIANAMASY

    double urand[PROZ];
    double urandalfa[PROZ];
    double Alpha22[PROZ];

    double flux2[PROZ], Teff2[PROZ], flux3[PROZ], Teff3[PROZ];
    double fluxcor[PROZ];
    double lum2, lum3, lumcor;

    double lasttimedynprom[PROZ];

    double curdt;

    RANDOMDEF randd[PROZ];

}
SOLINTIME;

typedef struct {
    double dsdt[PROZ];
    double dtempdt[PROZ];
    // WYMIANAMASY
    double dsdt_cor[PROZ];
    // !WYMIANAMASY
}
DERIVAT;

#define MAGROZ 500000
//#define MAGROZ 10000000

typedef struct {
    int    oczko;
    double rmag;
    double taudynmag;

    double lasttimedyn;
    double urand;
    double P1;
}
TABMAGSTRUCT;


typedef struct {
	int first;
	int last;
}TABMARANGE;

typedef struct {
    FORONERAD  tab[PROZ];
    TABMAGSTRUCT tabmag[MAGROZ];
//    TABMARANGE  tabmagRange[PROZ];
    TABMAGSTRUCT tabmagalfa[MAGROZ];
    TABMARANGE  tabmagalfaRange[PROZ];
    double dr[PROZ];
    double y[PROZ];
    double dy[PROZ];
    int IMMAX;
    int IMMAXalfa;
    double z[ZROZ];
}
GRID;



void derivatives(double Dt, SOLINTIME *curtm, DERIVAT *ddt);
void calculations(double Dt, SOLINTIME *curtm, SOLINTIME *newtm, DERIVAT *ddt, double derivscale, int calclum);

extern GRID gr;


extern double ptime;
extern int    ntime;

extern int rank;   //current thread no
extern int size;   //number of threads


#endif
