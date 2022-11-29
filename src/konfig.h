
#if !defined(KONFIG_H__INCLUDED_)
#define KONFIG_H__INCLUDED_


#define ln10   log(10)
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))


#define  Msun     ((double) 1.989e33)        /*constants in cgs units*/

#define  Hplanck  ((double) (6.6261e-27))
#define  Pi       ((double) (3.14159265358979323846))
#define  Sigmabol ((double) (5.6705e-5))
#define  Graw     ((double) (6.6726e-8))
#define  Crad     ((double) (2.9979e10))
#define  Mprot    ((double) (1.6726e-24))
#define  Year     ((double) ( 3.155693e7))
#define  Kbol     ((double) (1.3807e-16))


#define MROZ 100        /*number of points in stability curve*/
//#define PROZ 98        /*number of steps in radius*/
#define PROZ 218
//#define PROZ 386
//#define PROZ 49 /*for Fiamma model 4A */
//#define PROZ 196 /*for Fiamma model 4B */
#define ZROZ 99        /*size of z demiensions, for hdf only*/
//#define PROZ 402

typedef struct {

    double  Kappa;   // = 0.34;        /*parameters*/
    double  Alpha;   // = 1.e-2;
    double  AlphaCor;   // = 1.e-2;
    double  Mdot;   // = 1.0e-8;
    double  Qadv;   // = 0.33;

    double  BetaS;   // = 0.1;
    double  Kappa_dyn;   // = 100;

    double  Xicor;   // = 0.5;
    double  Dzet;   // = 0.0;
    double  etajet3; // = 1.e-7;


    double  Rmax;   // = 300.0;
    double  Rmin;   // = 3.01;

    double  Tmin;   // = 0.0;
    double  Tmax;   // = 1.e15;

    double  DtStart;   // = 1.e-4;
    int     krok_staly;   //0
    double  ErrMax;   // = 1;


    double  Mass;   // = (10*Msun);


    double  fac;   // = 0.02125;

    double  C1;   // = 1.25;       /*coefficients for vertical structure integration*/
    double  C2;   // = 6.25;
    double  C3;   // = 0.17;
    double  C4;   // = 1.0;
    double  C5;   // = 1.59;                /*for Xicor=0.5*/
    double  C6;   // = 1.0;

    double  B1;   // = 5.e-1;
    double  B2;   // = 5.e-1;

    double  xi0;   // = 8.0;              /*coefficient in modified viscosity law*/
    int     heating;   // = 0;            /*if 1 then corona in heating term -- if 0 then corona in cooling term*/
    int     viscosity;   // = 0;          /*if 1 then modified viscosity law*/
    int     hotvis;   // = 0;             /*if 1 then alpha hot < alpha cold*/
    double  Alpha0;   // = 0.8e-3;
    int     iptot;   // = 0;              /*if 1 then heating propto total pressure*/
    int     ipsqrt;   // = 1;              /*if 1 then heating propto sqrt(gas*rad) pressure*/
    int     PBK;   // = 1;                /*if 1 then first PBK_time sec without corona*/
    double  PBK_time;  //5.e5             /*time to activate corona*/

    int     wymiana_masy; // = 1;
    int 		wymiana_flux;
    char 		case_mz;
    int 		grid_radius;
    int gridHR;
    int gridRR;

    int     nsolve_method; //= 0     /*0 - Euler, 1 - RungeKutta, 2 - PredictorCorrector*/

    int     lekran;            //1000  /*wypisywanie na ekran*/

    int     llum;            //
    double  tlum;     // = 0.0   //Czas od ktorego beda zapisywane jasnosci
    double  dtlum;   // 10.0    czestotliwosc zapisywanie jasnosci
    double  deltalum;     // zmiana jasnosci powodujaca zapis do pliku

    int     lrad;            //
    double  trad;     // = 4.e5   //Czas od ktorego beda zapisywane wydruki radialne
    double  tradend;     // 4.5e5   //Czas do ktorego beda zapisywane wydruki radialne
    double  dtrad;   // 10.0    czestotliwosc zapisywanie wydrukow radialnych

    int     lsig;            //
    double  tsig;     // = 4
    double  dtsig;   // 10.0    czestotliwosc zapisywanie wydrukow radialnych

    int     lhdf;            //
    double  thdf;     // = 4
    double  dthdf;   // 10.0    czestotliwosc zapisywanie wydrukow radialnych
    int     lasthdf;


    int     lrestore;     //czy zapisywac dumpy do odtworzenia
    double  trestore;     // czas od ktorego beda zapisywane restory
    double  dtrestore;   //czestotliwosc zapisywania dumpor
    int     nrestore;    //> 0 muer dumpu do odtworzenia

    int     ji, jo;        //poczatek i koniec okna obliczniowego

    int     external_mdot;  //type of outer boundary
    int     alfa_flickering;  //changing viscosity with time and radius
    int       flicker_dyn;      //timescale of flickering scaled with dynamical time
    int       flicker_visc;     //timescale of flickering scaled with viscous time
    double bpar;             //scaling factor for viscosity flickering
    double  Tper;           //Tperturbacji

    char    rezerwa[100];  //for future use

} KONFIG;

int   InitDefaultKonfig(void);
void  ShowGreeting(void);       //Ekran powitalny
int   ReadConfig(int restart);



extern KONFIG  konf;

extern double  Mext;   // (Mdot*Msun/Year)
extern double  Rschw;  //  (2.0*Graw*Mass/Crad/Crad)
extern double  Ledd;   //  (1.25e39*(Mass/10.0/Msun))


#endif
