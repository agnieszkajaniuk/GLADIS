
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "konfig.h"
#include "dysk_zwm.h"
#include "konfig.h"
#include "readini.h"


KONFIG konf;

double  Mext;   // (Mdot*Msun/Year)
double  Rschw;  //  (2.0*Graw*Mass/Crad/Crad)
double  Ledd;   //  (1.25e39*(Mass/10.0/Msun))


/**************************************************************/
int InitDefaultKonfig(void)
/**************************************************************/
{
 /*parameters*/

    konf.Kappa = 0.34;       // Electron scattering opacity
    konf.Alpha = 1.e-2;      // Viscosity in the disk
    konf.AlphaCor = 1.e-2;   // Viscosity in the corona
    konf.Mdot = 1.0e-8;      // Accretion rate in Solar mass per year
    konf.Qadv = 0.33;        // Advection constant

    konf.BetaS = 0.1;		//Ratio of the disk magnetic field B_z to the global magnetic field (<=1.0)
    konf.Kappa_dyn = 100;	//Ratio of the magnetic dynamo timescale to the orbital timescale

    konf.Xicor = 0.5;		//Coronal fraction
    konf.Dzet = 0.0;		//Outflow strength
    konf.etajet3 = 0.0;		//Outflow fraction


    konf.Rmax = 300.0;        // Outer radius of the disk in Schwarzschild radii
    konf.Rmin = 3.01;         //Inner radius of the disk
    konf.Tmin = 0.0;          // Start time of evolution 
    konf.Tmax = 1.e15;        //End time of evolution, in seconds

    konf.DtStart = 1.e-4;     // Initial time-step
    konf.krok_staly = 0;      // If 1 then constant time step; if 0, then adaptive time-step

    konf.Mass = 10.0*Msun;    // Black hole mass in Solar masses


    konf.ErrMax = 1.0;

    konf.fac = 0.02125;

    konf.C1 = 1.25;       /*coefficients for vertical structure integration*/
    konf.C2 = 6.25;
    konf.C3 = 0.17;
    konf.C4 = 1.0;
    konf.C5 = 1.59;                /*for Xicor=0.5*/
    konf.C6 = 1.0;

    konf.B1 = 5.e-1;
    konf.B2 = 5.e-1;

    konf.xi0 = 8.0;              /*coefficient in modified viscosity law*/
    konf.heating = 0;            /*if 1 then corona in heating term -- if 0 then corona in cooling term*/
    konf.viscosity = 0;          /*if 1 then modified viscosity law*/
    konf.hotvis = 0;             /*if 1 then alpha hot < alpha cold*/
    konf.Alpha0 = 0.8e-3;
    konf.iptot = 1;                  /*if 1 then heating propto total pressure, if 0 then heating given by gas pressure. Works if ipsqrt is 0*/
    konf.ipsqrt = 0;                  /*if 1 then heating propto square root of total times gas pressure*/
    konf.PBK = 1;                /*if 1 then code initialized without corona*/
    konf.PBK_time = 5.e5;             /*time to activate corona, in seconds*/

    konf.wymiana_masy = 1;      /* if 1 then mass exchange between disk and corona*/
    konf.wymiana_flux = 0;      /*if 1 then mass evaporation is due to the generated flux. If 0 then mass evaporation due to magnetic field*/
    konf.case_mz = 'C';		/* Works if mass evaporation due to magnetic field. Options to choose: A - alvfen velocity,  B - light speed, C - alvfen devided by inner boundary condition */

    konf.grid_radius = 0;       /*if 1 then magnetic grid proportional to  square root of the disk radius. If 0 then two options to choose are as below*/

    konf.gridHR = 1;		/*magnetic cell proportional to disk radius*/
    konf.gridRR = 0;		/*magnetic cell proportional to disk height. Must be read from file rsig.dat*/

    konf.alfa_flickering = 0;	/* If 1 - changing viscosity with time and radius as a Markov chain. Caution: file rsig.dat needed! ;strength of this flickering governed by the bpar parameter */

    konf.flicker_dyn = 1;	/*If 1 - timescale of flickering scaled with dynamical time */
    konf.flicker_visc = 0;	/*If 1 - timescale of flickering scaled with viscous time */

    konf.bpar = 0.1;	/* scaling factor for the viscosity flickering */

    konf.nsolve_method = 0; //Numerical method choice    /*0 - Euler, 1 - RungeKutta, 2 - PredictorCorrector, 3 - Heun*/

    konf.lekran = 1000;      // Screen printout frequency, in timesteps

    konf.llum = 1;        // If 1 then save the lightcurve in lum.dat
    konf.tlum = 0.0;     //    Starting time for the lightcurve dumping 
    konf.dtlum = 10.0;   //    dumping frequency of the lightcurve
    konf.deltalum = 0.001;     // minimum relative change in the luminosity for which the luminisity is dumped in output file. If 0.0 then no extra output

    konf.lrad = 0;         // If 1 then save the radial  profiles to files rsig.*
    konf.trad = 5.0e4;     // Starting time for the radial profiles dumping
    konf.tradend = 6.0e4;     // End time for the radial profiles dumping
    konf.dtrad = 1000.0;   // Dumping frequency for the radial profiles

    konf.lsig = 0;         // If 1 then save the S-curves to files sigte.*
    konf.tsig = 0.0;        // Starting time for the S-curves dumping
    konf.dtsig = 10.0;   // Dumping frequency for the S-curves

    konf.lhdf = 0;         // If 1 then save the output in hdf format to files hdf*.*
    konf.thdf = 0.0;       // Starting time for the hdf files dumping
    konf.dthdf = 10.0;   // Dumping frequency for the hdf files
    konf.lasthdf = 0;    // End time for the hdf files dumping

    konf.lrestore = 0;     //If 1 then save restart files
    konf.trestore = 0.0;     // Starting time for the restart files dumping
    konf.dtrestore = 3600.0;   //Dumping frequency for the restart files (realtime seconds)
    konf.nrestore = 0;    // Number of restart file to be used, if > 0 then restoring evolution from dump.dat.nrestore

    konf.alfa_flickering = 0; // If 1 - changing viscosity with time and radius as a Markov chain. Caution: file rsig.dat needed! ;strength of this flickering governed by the bpar parameter
    konf.flicker_dyn = 1;	//If 1 - timescale of flickering scaled with dynamical time
    konf.flicker_visc = 0;	// If 1 - timescale of flickering scaled with viscous time
    konf.bpar = 0.1;		// Scaling factor for the flickering of alfa
    konf.external_mdot = 0;	// type of outer boundary. If 0 - constant mdot, 1 - periodic mdot, 2 - tabularized mdot
    konf.Tper = 0.0;		// Perturbation period if periodic outer boundary


    Mext   = konf.Mdot*Msun/Year;
    Rschw  = 2.0*Graw*konf.Mass/Crad/Crad;
    Ledd   = 1.25e39*(konf.Mass/10.0/Msun);


    return 1;
}



/**************************************************************/
void ShowGreeting(void)
//Ekran powitalny
/**************************************************************/
{

  

    printf("\n");
    printf("\tGGGG    L         A     DDDD   I   SSSS\n");
    printf("\tG       L        A A    D   D  I  S   \n");
    printf("\tG  GG   L       A   A   D   D  I   SSS\n");
    printf("\tG   G   L      AAAAAAA  D   D  I      S\n");
    printf("\tGGGG    LLLL  A       A DDDD   I  SSSS\n");

    printf("\n\n");
    printf("\t\tDEVELOPED BY AGNIESZKA JANIUK\n");
    printf("\t\tGLADIS V2.9.4 - 15/12/20\n");
    printf("\n\n");
    printf("\tCONFIGURATION:\n\n");


    printf("\tLedd=%E, Medd=%E\n", Ledd, Ledd/Crad/Crad*16.0);
    printf("\tMass=%E, mdot=%E\n", konf.Mass, Mext/(Ledd/Crad/Crad*16.0));
    printf("\tLgener=%E\n", Mext*Crad*Crad/16);
    printf("\tAlpha=%E, AlphaCor=%E\n", konf.Alpha, konf.AlphaCor);
    printf("\talfa_flickering=%d, bpar=%E\n", konf.alfa_flickering, konf.bpar);
    printf("\tflicker_dyn=%d\n", konf.flicker_dyn);
    printf("\tXicor=%E, Dzet=%E\n", konf.Xicor, konf.Dzet);
    printf("\tRschw=%E\n", Rschw);


    //     getc(stdin);


    if(!konf.ipsqrt) {
        if(!konf.iptot)
            printf("\tHeating proportional to Pgas\n");
        else
            printf("\tHeating proportional to Ptot\n");
    } else {
        printf("\tHeating proportional to sqrt(Pgas*Ptot)\n");
    }

    if(konf.wymiana_masy) {
        if(konf.PBK)
            printf("\tPBK: Mass exchange with corona:\n");
        printf("\tcorona after disk relaxation time t=%E\n", konf.PBK_time);

        if(!konf.wymiana_flux)
            printf("\tMass evaporation: magnetic dynamo\n");
        else
            printf("\tMass evaporation: flux generated in the disk\n");

        if(konf.case_mz == 'A')
            printf("\tAccretion rate in vertical direction: Alvfen velocity\n");
        if(konf.case_mz == 'B')
            printf("\tAccretion rate in vertical direction: Speed of light\n");
        if(konf.case_mz == 'C')
            printf("\tAccretion rate in vertical direction: Alvfen velocity devided by inner boundary condition\n");

        printf("\tBetaS = %E, Kappa_dyn = %E\n", konf.BetaS, konf.Kappa_dyn);
        //   getc(stdin);
    }

    printf("\tNumber of threads: %d\n", size);

    printf("\n\tCompiler dependencies:\n");
    printf("\tInteger size - %lu\n", sizeof(int));
    printf("\tLong size    - %lu\n", sizeof(long));
    printf("\tDouble size  - %lu\n", sizeof(double));


    printf("\n\n\n");

    //    getc(stdin);


}


/**************************************************************/
int ReadConfig(int restart)
/**************************************************************/
{
    const char *inifile = "dysk_zwm.ini";
    char *ptr = NULL;


    //gen
    if(!restart) {
        ptr = PathIni(inifile, "gen", "Mass");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Mass");
            return 0;
        }
        konf.Mass = Msun*strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Kappa");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Kappa");
            return 0;
        }
        konf.Kappa = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Alpha");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Alpha");
            return 0;
        }
        konf.Alpha = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "AlphaCor");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing%s\n", inifile, "AlphaCor");
            return 0;
        }
        konf.AlphaCor = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Mdot");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing%s\n", inifile, "Mdot");
            return 0;
        }
        konf.Mdot = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Qadv");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Qadv");
            return 0;
        }
        konf.Qadv = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "BetaS");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "BetaS");
            return 0;
        }
        konf.BetaS = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Kappa_dyn");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Kappa_dyn");
            return 0;
        }
        konf.Kappa_dyn = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Xicor");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Xicor");
            return 0;
        }
        konf.Xicor = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Dzet");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Dzet");
            return 0;
        }
        konf.Dzet = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "etajet3");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "etajet3");
            return 0;
        }
        konf.etajet3 = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "fac");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "fac");
            return 0;
        }
        konf.fac = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "C1");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "C1");
            return 0;
        }
        konf.C1 = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "C2");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "C2");
            return 0;
        }
        konf.C2 = strtod(ptr,NULL);

        ptr = PathIni(inifile, "gen", "C3");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "C3");
            return 0;
        }
        konf.C3 = strtod(ptr,NULL);

        ptr = PathIni(inifile, "gen", "C4");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "C4");
            return 0;
        }
        konf.C4 = strtod(ptr,NULL);

        ptr = PathIni(inifile, "gen", "C5");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "C5");
            return 0;
        }
        konf.C5 = strtod(ptr,NULL);

        ptr = PathIni(inifile, "gen", "C6");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "C6");
            return 0;
        }
        konf.C6 = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "B1");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "B1");
            return 0;
        }
        konf.B1 = strtod(ptr,NULL);

        ptr = PathIni(inifile, "gen", "B2");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "B2");
            return 0;
        }
        konf.B2 = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "xi0");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "xi0");
            return 0;
        }
        konf.xi0 = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Alpha0");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Alpha0");
            return 0;
        }
        konf.Alpha0 = strtod(ptr,NULL);


        ptr = PathIni(inifile, "gen", "Tper");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Tper");
            return 0;
        }
        konf.Tper = strtod(ptr,NULL);


        //switches

        ptr = PathIni(inifile, "switches", "heating");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "heating");
            return 0;
        }
        konf.heating = atoi(ptr);


        ptr = PathIni(inifile, "switches", "viscosity");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "viscosity");
            return 0;
        }
        konf.viscosity = atoi(ptr);

        ptr = PathIni(inifile, "switches", "hotvis");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "hotvis");
            return 0;
        }
        konf.hotvis = atoi(ptr);

        ptr = PathIni(inifile, "switches", "iptot");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "iptot");
            return 0;
        }
        konf.iptot = atoi(ptr);

        ptr = PathIni(inifile, "switches", "ipsqrt");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "ipsqrt");
            return 0;
        }
        konf.ipsqrt = atoi(ptr);

        ptr = PathIni(inifile, "switches", "PBK");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "PBK");
            return 0;
        }
        konf.PBK = atoi(ptr);

        ptr = PathIni(inifile, "switches", "PBK_time");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "PBK_time");
            return 0;
        }
        konf.PBK_time = strtod(ptr,NULL);




        ptr = PathIni(inifile, "switches", "wymiana_masy");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "wymiana_masy");
            return 0;
        }
        konf.wymiana_masy = atoi(ptr);



        ptr = PathIni(inifile, "switches", "wymiana_flux");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "wymiana_flux");
            return 0;
        }
        konf.wymiana_flux = atoi(ptr);


        ptr = PathIni(inifile, "switches", "case_mz");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "case_mz");
            return 0;
        }
        konf.case_mz = *ptr;

        ptr = PathIni(inifile, "switches", "grid_radius");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "grid_radius");
            return 0;
        }
        konf.grid_radius = atoi(ptr);


        ptr = PathIni(inifile, "switches", "gridRR");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "gridRR");
            return 0;
        }
        konf.gridRR = atoi(ptr);

        ptr = PathIni(inifile, "switches", "gridHR");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "gridHR");
            return 0;
        }
        konf.gridHR = atoi(ptr);

        if(konf.grid_radius == 0 && konf.gridRR == 0 && konf.gridHR == 0) {
            printf("Wrong file %s. Choose an option for magnetic grid %s %s %s \n", inifile, "grid_radius", "gridRR","gridHR" );
        }

        ptr = PathIni(inifile, "switches", "alfa_flickering");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "alfa_flickering");
            return 0;
        }

        konf.alfa_flickering = atoi(ptr);


        ptr = PathIni(inifile, "switches", "flicker_dyn");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "flicker_dyn");
            return 0;
        }

        konf.flicker_dyn = atoi(ptr);

        ptr = PathIni(inifile, "switches", "flicker_visc");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "flicker_visc");
            return 0;
        }

        konf.flicker_visc = atoi(ptr);


        ptr = PathIni(inifile, "switches", "bpar");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "bpar");
            return 0;
        }

        konf.bpar = atof(ptr);

        ptr = PathIni(inifile, "switches", "external_mdot");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "external_mdot");
            return 0;
        }

        konf.external_mdot = atoi(ptr);

        //grid

        ptr = PathIni(inifile, "grid", "Rmax");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Rmax");
            return 0;
        }
        konf.Rmax = strtod(ptr,NULL);

        ptr = PathIni(inifile, "grid", "Rmin");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Rmin");
            return 0;
        }
        konf.Rmin = strtod(ptr,NULL);


        //solve

        ptr = PathIni(inifile, "solve", "Tmin");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "Tmin");
            return 0;
        }
        konf.Tmin = strtod(ptr,NULL);

    } //if(!restart)

    ptr = PathIni(inifile, "solve", "Tmax");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "Tmax");
        return 0;
    }
    konf.Tmax = strtod(ptr,NULL);

    if(!restart) {

        ptr = PathIni(inifile, "solve", "DtStart");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "DtStart");
            return 0;
        }
        konf.DtStart = strtod(ptr,NULL);


        ptr = PathIni(inifile, "solve", "krok_staly");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "krok_staly");
            return 0;
        }
        konf.krok_staly = atoi(ptr);

    } //if(!restart)


    ptr = PathIni(inifile, "solve", "ErrMax");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "ErrMax");
        return 0;
    }
    konf.ErrMax = strtod(ptr,NULL);


    if(!restart) {
        ptr = PathIni(inifile, "solve", "nsolve_method");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "nsolve_method");
            return 0;
        }
        konf.nsolve_method = atoi(ptr);
    } //if(!restart)

    //inout

    ptr = PathIni(inifile, "inout", "lekran");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "lekran");
        return 0;
    }
    konf.lekran = atoi(ptr);

    if(!restart) {

        ptr = PathIni(inifile, "inout", "llum");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "llum");
            return 0;
        }
        konf.llum = atoi(ptr);

        ptr = PathIni(inifile, "inout", "tlum");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "tlum");
            return 0;
        }
        konf.tlum = strtod(ptr,NULL);

    } //if(!restart)

    ptr = PathIni(inifile, "inout", "dtlum");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "dtlum");
        return 0;
    }
    konf.dtlum = strtod(ptr,NULL);

    ptr = PathIni(inifile, "inout", "deltalum");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "deltalum");
        return 0;
    }
    konf.deltalum = strtod(ptr,NULL);

    if(!restart) {

        if (konf.llum > 0)
            konf.tlum -= konf.dtlum;


        ptr = PathIni(inifile, "inout", "lrad");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "lrad");
            return 0;
        }
        konf.lrad = atoi(ptr);

        ptr = PathIni(inifile, "inout", "trad");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "trad");
            return 0;
        }
        konf.trad = strtod(ptr,NULL);

        ptr = PathIni(inifile, "inout", "tradend");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "tradend");
            return 0;
        }
        konf.tradend = strtod(ptr,NULL);

    } //    if(!restart)



    ptr = PathIni(inifile, "inout", "dtrad");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "dtrad");
        return 0;
    }
    konf.dtrad = strtod(ptr,NULL);


    if(!restart) {

        if (konf.lrad > 0)
            konf.trad -= konf.dtrad;

        ptr = PathIni(inifile, "inout", "lsig");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "lsig");
            return 0;
        }
        konf.lsig = atoi(ptr);

        ptr = PathIni(inifile, "inout", "tsig");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "tsig");
            return 0;
        }
        konf.tsig = strtod(ptr,NULL);
    } //    if(!restart)


    ptr = PathIni(inifile, "inout", "dtsig");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "dtsig");
        return 0;
    }
    konf.dtsig = strtod(ptr,NULL);


    if(!restart) {

        if (konf.lsig > 0)
            konf.tsig -= konf.dtsig;

        ptr = PathIni(inifile, "inout", "lhdf");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "lhdf");
            return 0;
        }
        konf.lhdf = atoi(ptr);

        ptr = PathIni(inifile, "inout", "thdf");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "thdf");
            return 0;
        }
        konf.thdf = strtod(ptr,NULL);
    } //    if(!restart)


    ptr = PathIni(inifile, "inout", "dthdf");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "dthdf");
        return 0;
    }
    konf.dthdf = strtod(ptr,NULL);


    if (!restart) {

        if (konf.lhdf > 0)
            konf.thdf -= konf.dthdf;
        konf.lasthdf = 0;

        //restore
        ptr = PathIni(inifile, "restore", "lrestore");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "lrestore");
            return 0;
        }
        konf.lrestore = atoi(ptr);


        ptr = PathIni(inifile, "restore", "trestore");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "trestore");
            return 0;
        }
        konf.trestore = strtod(ptr,NULL);

    } //if(!restart)

    ptr = PathIni(inifile, "restore", "dtrestore");
    if (ptr == NULL) {
        printf("Wrong file %s. Parameter missing %s\n", inifile, "dtrestore");
        return 0;
    }
    konf.dtrestore = strtod(ptr,NULL);

    if(!restart) {

        if (konf.lrestore > 0)
            konf.trestore += time(NULL) - konf.dtrestore;

        ptr = PathIni(inifile, "restore", "nrestore");
        if (ptr == NULL) {
            printf("Wrong file %s. Parameter missing %s\n", inifile, "nrestore");
            return 0;
        }
        konf.nrestore = atoi(ptr);


        Mext   = konf.Mdot*Msun/Year;
        Rschw  = 2.0*Graw*konf.Mass/Crad/Crad;
        Ledd   = 1.25e39*(konf.Mass/10.0/Msun);
    } //    if(!restart)

    return 1;
}










































