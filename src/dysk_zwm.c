#ifdef MPI_USED
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <fcntl.h>


#include "konfig.h"
#include "dysk_zwm.h"
#include "nsolve.h"
#include "inout.h"
#include "steady.h"
#include "random.h"
#include "boundary.h"


GRID gr;

double ptime;
int    ntime;

int rank = 0;   //current thread no
int size = 1;   //number of threads


/*functions*/
void    InitializeDisk(SOLINTIME *curtm);
void    derivatives(double Dt, SOLINTIME *curtm, DERIVAT *ddt);
void    calculations(double Dt, SOLINTIME *curtm, SOLINTIME *newtm, DERIVAT *ddt, double derivscale, int calclum);
void    ProgramStop(int code);

double 	Markoff(double urand, RANDOMDEF *rd);
static inline void 		bval1(SOLINTIME *curtm);
static inline void    bval2(SOLINTIME *curtm);

double etajet;

/**************************************************************/
int main (int argc, char **argv)
/**************************************************************/
{

    double dt;

    int    res;
    SOLINTIME curtm;      //current solution
    SOLINTIME newtm;      //new solution


    memset(&curtm, 0, sizeof(SOLINTIME));
    memset(&newtm, 0, sizeof(SOLINTIME));

    //Obsluga sysgnalow

    //    signal(SIGINT, ProgramStop);   //Ctrl-C
    //    signal(SIGQUIT, SIG_IGN);
    //    signal(SIGHUP, SIG_IGN);
    //    signal(SIGTERM, SIG_IGN);

    //    signal(SIGTERM, ProgramStop);     //kill
    //    signal(SIGCHLD, SIG_IGN);



    //Ustawienie zerowej dlugosci buforow std
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

#ifndef WIN32
    fcntl(fileno(stdout),F_SETFL,O_SYNC);
    fcntl(fileno(stderr),F_SETFL,O_SYNC);
#endif
    InitDefaultKonfig();  //Domyslna konfiguracja

    if(!ReadConfig(0)) {   ////Wczytanie konfiguracji
        ProgramStop(0);
        return 1;
    }


#ifdef MPI_USED

    int ierr;
    ierr = MPI_Init(&argc, &argv);          //Start MPI


    MPI_Comm_rank( MPI_COMM_WORLD, &rank );  //Pobranie numeru procesu, 0 - proces glowny
    MPI_Comm_size( MPI_COMM_WORLD, &size );  //Pobranie liczby procesow

    if(((PROZ-2)%size) != 0) {
        printf("Number of threads and PROZ dosn't match!\n");
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    int slabsize = (PROZ-2)/size;

    konf.ji = 1 + slabsize*rank;
    konf.jo = slabsize*(rank+1);

    if(size > 1 && (konf.viscosity || konf.hotvis))
    {
        printf("Warning: models with modified viscosity!\n");
	//        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

#else

    konf.ji = 1;
    konf.jo = PROZ - 2;

#endif

    if(rank == 0)
        ShowGreeting();       //Ekran powitalny


    if(konf.nrestore > 0) { //odczytanie dumpa

        printf("Restarting from dump nr %d rank %d...\n", konf.nrestore, rank);
        if(RestartRestore(konf.nrestore, &curtm, &dt)) {
            printf("Dump nr %d rank %d restarted \n", konf.nrestore, rank);
        } else {

            if(konf.nrestore == 1) {
                printf("Restarting problem...\n");
                konf.nrestore = 0;
                goto init;
            }

            printf("Cannot restart from dump nr %d rank %d.\n", konf.nrestore, rank);
            return 1;
        }

        //Ponowne odcztanie konfiguracji
        if(!ReadConfig(1)) {   ////Doczytanie z konfiguracji zmienionych wart. (tylko delty zapisu do pliku)
            ProgramStop(0);
            return 1;
        }

        if (konf.lrestore > 0)
            konf.trestore = time(NULL);

    } else {  //inicjalizacja dysku
init:
        InitializeDisk(&curtm);

        ntime = 0;
        ptime = konf.Tmin;
        dt = konf.DtStart;

        WriteKrzyweS();
    }

    if(!InitializeOutput(konf.nrestore)) {
        ProgramStop(0);
        return 1;
    }

    if(konf.nrestore == 0) { //Pierwszy dump restartu zaraz po obliczniu modelu stacjonarnego
        if (konf.lrestore > 0 && time(NULL) >= konf.trestore + konf.dtrestore) {
            konf.trestore += konf.dtrestore;
            RestartDump(&curtm, &dt);
        }
    }

    if(rank == size - 1)
        InitOuterBoundary();

    if(rank == 0) {
        printf("Evolution begins...\n");
        printf("Time = %lf\n", ptime);
    }

	clock_t lastClock = clock();
  
    
    while (ptime <= konf.Tmax)   /*time evolution*/
    {

        newtm.curdt = dt;

//        Boundary(&curtm);           //Apply boundary conditions, tests and review are needed

        if(rank == size - 1)
            OuterBoundary(&Mext, &curtm, ptime);           //Apply outer boundary conditions


        switch (konf.nsolve_method) {

        case 0:
            res = Euler(&dt, &curtm, &newtm);
            break;
        case 1:
            res = RungeKutta(&dt, &curtm, &newtm);
            break;
        case 2:
            /*Adams-Moulton Predictor-Corrector Method with variable step-size*/
            res = PredictorCorrector(&dt, ntime, &curtm, &newtm);
            break;
        case 3:
            res = Heun(&dt, &curtm, &newtm);
            break;
        default:
            res = 0;
            if(rank == 0)
                printf("Method undefined!\n");
            ProgramStop(0);
            break;
        }


        if(!res)
            continue;

        ptime += dt;
        ntime++;

        if((ntime % konf.lekran) == 0) {
            double lum2 = newtm.lum2;
            double lumcor = newtm.lumcor;

#ifdef MPI_USED
            //Pobieramy sume jasnosci
            if(size > 1) {
                double isum[2], osum[2];
                isum[0] = lum2;
                isum[1] = lumcor;

                MPI_Reduce( isum, osum, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

                lum2   = osum[0];
                lumcor = osum[1];
            }

            MPI_Bcast ( &Mext, 1, MPI_DOUBLE, size - 1, MPI_COMM_WORLD );
#endif
            if(rank == 0) {

                lum2 = log10(lum2);

                if(konf.wymiana_masy) {
                    lumcor = log10(lumcor);
                    printf("time=%E, Dt= %E, lum=%E lumcor=%E,  Mext=%E calcTime=%.3lf sec\n", ptime, dt, lum2, lumcor, Mext, ((double)(clock()-lastClock))/CLOCKS_PER_SEC);
                } else
                    printf("time=%E, Dt= %E, lum=%E, etajet=%E, Mext=%E calcTime=%.3lf sec\n", ptime, dt, lum2, etajet, Mext, ((double)(clock()-lastClock))/CLOCKS_PER_SEC);
			    lastClock = clock();                    
            }
        }

        dataio(ptime, ntime, dt, &newtm);

        memcpy(&curtm, &newtm, sizeof(SOLINTIME));

    }

    if(rank == 0)
        printf("Evolution finished\n, Time=%E\n", ptime);

    ProgramStop(0);
}

/**************************************************************/
void InitializeDisk(SOLINTIME *curtm)
/**************************************************************/
{

    register int jr;
    register int i,j=0;
    FILE *fcachestat;
    char s_file[256];
    double r1, h1, s1, s2, s3, s4, s5, s6;
    double ksi;

    //siatka
    if(rank == 0)
        printf("Initializing disk: grid y; %d points\n", PROZ);

    gr.tab[0].r = konf.Rmin;                         //siatka promieni
    gr.tab[0].rr = konf.Rmin;
    gr.tab[PROZ-1].r =konf.Rmax;

    gr.dr[0] = konf.fac;
    gr.dr[PROZ-1] = (gr.tab[PROZ-1].r - gr.tab[0].r)*konf.fac+konf.fac;

    gr.y[0]=2*sqrt(konf.Rmin*Rschw);                      //pomocnicza siatka do rownan ewolucji
    gr.y[PROZ-1]=2*sqrt(konf.Rmax*Rschw);                 //promien w cm!!

    for (i=1; i<PROZ; i++) {

        gr.dy[i]=(gr.y[PROZ-1]-gr.y[0])/(PROZ-1);
        gr.y[i]=gr.y[0]+gr.dy[i]*i;
        gr.tab[i].r=pow((gr.y[i]/2.0), 2)/Rschw;
        gr.dr[i]=gr.tab[i].r-gr.tab[i-1].r;
    }

    gr.dy[0]=gr.dy[1];


    //  printf(" r[0]:  r=%E, dr=%E, y=%E, dy=%E\n", gr.tab[0].r, gr.dr[0], gr.y[0], gr.dy[0]);

	//Grud Z
	{
		double	maxHd = konf.Rmax*Rschw/10.0;
		double	ratio = 1.44;
		int		nzones = ZROZ/2;

		gr.z[ZROZ/2] = 0.0;
		
		double dxfac = (maxHd - 0)*(ratio - 1.0)/(pow(ratio, nzones) - 1.0)*ratio;
		for ( i = ZROZ/2+1; i < ZROZ; i++ ){
				gr.z[i] = gr.z[i-1] + dxfac;
				dxfac *= ratio;
		}

		dxfac = (maxHd - 0)*(ratio - 1.0)/(pow(ratio, nzones) - 1.0)*ratio;
		for ( i = ZROZ/2-1; i >= 0; i-- ){
				gr.z[i] = gr.z[i+1] - dxfac;
				dxfac *= ratio;
		}
	}
/*	
    for ( i = 0; i < ZROZ; i++ ){
		printf("i=%03d z=%E maxHd=%E\n", i, gr.z[i], maxHd);
	}
*/

	
	
    ///////////////////////////////////////
    /* Initializing grid for magnetic dynamo*/

    if(konf.grid_radius)
        if(rank == 0)
            printf("Magnetic grid: dRmag = Sqrt(R) \n");

    //magnetic grid dr = H
    if(!konf.grid_radius) {

        if(konf.gridHR) {

            sprintf(s_file, "rsig.dat");
            if(rank == 0)
                printf("Reading file rsig.dat \n");

            if((fcachestat = fopen(s_file,"r"))==NULL) {
                if(rank == 0)
                    printf("Cannot open rsig.dat\n");
                exit(0);
            }

            if(rank == 0)
                printf("Magnetic grid: dRmag = H \n");

            //   gr.tab[0].Rmag = gr.tab[0].r*Rschw;

            for (i=0; i<PROZ; i++) {
                if(fscanf(fcachestat,"%lf %lf %lf %lf %lf\n", &r1, &s1, &s2, &s3, &h1)!= 5) {
                    //  if(rank == 0)
                    printf("fscanf eror 1\n");
                    getc(stdin);
                }


                //    gr.tab[i].Rmag = gr.tab[i-1].Rmag+h1;
                gr.tab[i].Rmag = h1*0.98;

                //            printf("%E %E\n", gr.tab[i].Rmag, h1);
                //            printf("Rmag=%E cm, %E Rschw, dR=%E\n", gr.tab[i].Rmag, gr.tab[i].Rmag/Rschw);

            }
            fclose(fcachestat);
        }


        if(konf.gridRR) {

            if(rank == 0)
                printf("Magnetic grid: dRmag = 0.05 R \n");

            for (i=1; i<PROZ-1; i++) {
                gr.tab[i].Rmag = gr.tab[i].r*Rschw*0.05;
                //	printf("r=%E dR=%E\n", gr.tab[i].r, gr.tab[i].Rmag );
            }
        }

        gr.tab[0].Rmag = -1.0*(gr.tab[1].Rmag*(gr.tab[0].r*Rschw-gr.tab[1].r*Rschw)
                               +gr.tab[2].Rmag*(gr.tab[1].r*Rschw - gr.tab[0].r*Rschw))/(gr.tab[1].r*Rschw - gr.tab[2].r*Rschw);
        gr.tab[PROZ-1].Rmag = -1.0*(gr.tab[PROZ-2].Rmag*(gr.tab[PROZ-1].r*Rschw-gr.tab[PROZ-3].r*Rschw)
                                    -gr.tab[PROZ-3].Rmag*(gr.tab[PROZ-1].r*Rschw - gr.tab[PROZ-2].r*Rschw))
                              /(gr.tab[PROZ-2].r*Rschw - gr.tab[PROZ-1].r*Rschw);


        if(rank == 0)
            printf(" Rmag_0 = %E Rmag_PROZ = %E \n", gr.tab[0].Rmag, gr.tab[PROZ-1].Rmag);


        int im, ii, jj;
        double hmag;


        im = 0;
        gr.tabmag[0].rmag = gr.tab[0].r*Rschw;
        gr.tabmag[0].oczko = 0;


        //       printf("%E, %E\n", gr.tabmag[im].rmag, gr.tab[PROZ-1].Rmag);
        //			printf("Rschw=%E\n",Rschw);

        while (gr.tabmag[im].rmag < gr.tab[PROZ-1].r*Rschw) {
            //    printf("Tutaj jestem 1, im=%d\n", im);

            for (jj = 0; jj <= PROZ-1 ; jj++) {
                if(gr.tabmag[im].rmag < gr.tab[jj+1].r*Rschw)
                    break;
            }

            //            printf("jj=%d, H1=%E, H2=%E r1=%E, r2=%E \n", jj, gr.tab[jj].Rmag, gr.tab[jj+1].Rmag, gr.tab[jj].r*Rschw, gr.tab[jj+1].r*Rschw );


            hmag = gr.tab[jj].Rmag + (gr.tab[jj+1].Rmag-gr.tab[jj].Rmag)*(gr.tabmag[im].rmag - gr.tab[jj].r*Rschw)/(gr.tab[jj+1].r*Rschw - gr.tab[jj].r*Rschw);

            //            	printf("im=%d, hmag = %E jj=%d, akt=%E [%E, %E]\n", im, hmag, jj, gr.tabmag[im].rmag, gr.tab[jj].r*Rschw, gr.tab[jj+1].r*Rschw);

            im++;
            gr.tabmag[im].rmag =   gr.tabmag[im-1].rmag + hmag;
            gr.tabmag[im].oczko = jj;

            //  printf("%E\n", gr.tabmag[im].rmag);
            //    getc(stdin);

        }


        gr.IMMAX = im;

        if(rank == 0)
            printf("gr.IMMAX = %d\n", gr.IMMAX);

        if(gr.IMMAX >= MAGROZ -1) {
            if(rank == 0)
                printf("gr.IMMAX > MAGROZ\n");
            ProgramStop(0);
            return;
        }


        //       getc(stdin);

        for (ii=0; ii<=gr.IMMAX; ii++) {

            double Omegamag = pow(Graw*konf.Mass/pow(gr.tabmag[ii].rmag,3.0), 0.5);
            gr.tabmag[ii].taudynmag = konf.Kappa_dyn/Omegamag;

            gr.tabmag[ii].lasttimedyn = 0.0;
            gr.tabmag[ii].urand = 0.0;
            //            printf("ii=%d, rmag=%E, taudynmag=%E \n", ii, gr.tabmag[ii].rmag,
            //                   gr.tabmag[ii].taudynmag);
        }


    } //if(!konf.grid_radius)

    ////////////////////////////////////////////
    /* Initializing grid for alpha flickering */
    if(konf.alfa_flickering) {
        {
            sprintf(s_file, "rsig.dat");
            if(rank == 0)
                printf("Reading file rsig.dat for alfa fickering\n");

            if((fcachestat = fopen(s_file,"r"))==NULL) {
                if(rank == 0)
                    printf("Cannot open rsig.dat\n");
                exit(0);
            }

            if(rank == 0)
                printf("Flickering alpha grid: dR = H \n");

            //   gr.tab[0].Rmag = gr.tab[0].r*Rschw;

            for (i=0; i<PROZ; i++) {
	      if(fscanf(fcachestat,"%lf %lf %lf %lf %lf %lf %lf %lf\n", &r1, &s1, &h1, &s2, &s3, &s4, &s5, &s6)!= 8) {
//a                    if(rank == 0)
                    printf("fscanf eror 2\n");
                    getc(stdin);
                }


                //    gr.tab[i].Rmag = gr.tab[i-1].Rmag+h1;

                printf("r=%lf, R=%E,  h=%E\n", r1, r1*Rschw, h1);
                gr.tab[i].Rmagalfa = h1*0.98;

            }

            fclose(fcachestat);


            for (i=0; i<PROZ-1; i++) {
                printf("i=%d, Rmagalfa=%E\n", i, gr.tab[i].Rmagalfa);
            }

            gr.tab[0].Rmagalfa = -1.0*(gr.tab[1].Rmagalfa*(gr.tab[0].r*Rschw-gr.tab[1].r*Rschw)
                                       +gr.tab[2].Rmagalfa*(gr.tab[1].r*Rschw - gr.tab[0].r*Rschw))/(gr.tab[1].r*Rschw - gr.tab[2].r*Rschw);


            gr.tab[PROZ-1].Rmagalfa = -1.0*(gr.tab[PROZ-2].Rmagalfa*(gr.tab[PROZ-1].r*Rschw-gr.tab[PROZ-3].r*Rschw)
                                            -gr.tab[PROZ-3].Rmagalfa*(gr.tab[PROZ-1].r*Rschw - gr.tab[PROZ-2].r*Rschw))
                                      /(gr.tab[PROZ-2].r*Rschw - gr.tab[PROZ-1].r*Rschw);


            if(rank == 0)
                printf(" Rmagalfa_0 = %E Rmagalfa_PROZ = %E \n", gr.tab[0].Rmagalfa, gr.tab[PROZ-1].Rmagalfa);


            int im, ii, jj;
            double hmag;


            im = 0;
            gr.tabmagalfa[0].rmag = gr.tab[0].r*Rschw;
            gr.tabmagalfa[1].rmag = gr.tab[0].r*Rschw;
            gr.tabmagalfa[0].oczko = 0;

            for (jj = 0; jj < PROZ ; jj++) {
            	gr.tabmagalfaRange[jj].first = -1;
            	gr.tabmagalfaRange[jj].last = -1;
            }
            gr.tabmagalfaRange[jj].first = 0;
            
            printf("Siatka %E - %E\n", gr.tab[0].r*Rschw, gr.tab[PROZ-1].r*Rschw);

            while (gr.tabmagalfa[im].rmag < gr.tab[PROZ-1].r*Rschw) {
                //	    printf("1, im=%d, rmag=%E \n", im, gr.tabmagalfa[im].rmag);

                for (jj = 0; jj < PROZ-1; jj++) {
                    if(gr.tabmagalfa[im].rmag < gr.tab[jj+1].r*Rschw)
                        break;
                }

                hmag = gr.tab[jj].Rmagalfa + (gr.tab[jj+1].Rmagalfa-gr.tab[jj].Rmagalfa)*(gr.tabmagalfa[im].rmag - gr.tab[jj].r*Rschw)/(gr.tab[jj+1].r*Rschw - gr.tab[jj].r*Rschw);

                //	    printf("jj=%d, Rmagalfa=%E, hmag = %E\n", jj, gr.tab[jj].Rmagalfa, hmag);
                //	    printf("delta = %E, skladniki=%E, %E, delta1=%E, delta2=%E, delta3=%E\n", (gr.tab[jj+1].Rmagalfa-gr.tab[jj].Rmagalfa)*(gr.tabmagalfa[im].rmag - gr.tab[jj].r*Rschw)/(gr.tab[jj+1].r*Rschw - gr.tab[jj].r*Rschw),

                //		   gr.tab[jj+1].Rmagalfa, gr.tab[jj].Rmagalfa,
                //	   gr.tab[jj+1].Rmagalfa-gr.tab[jj].Rmagalfa,
                //	   gr.tabmagalfa[im].rmag - gr.tab[jj].r*Rschw,
                //gr.tab[jj+1].r*Rschw - gr.tab[jj].r*Rschw);

                im++;
                gr.tabmagalfa[im].rmag =   gr.tabmagalfa[im-1].rmag + hmag;
                gr.tabmagalfa[im].oczko = jj;
                
                if(gr.tabmagalfaRange[jj].first == -1) gr.tabmagalfaRange[jj].first = im;
                if(gr.tabmagalfaRange[jj].last < im) gr.tabmagalfaRange[jj].last = im;
            }

//            for (jj = 1; jj < PROZ-1; jj++) {
//            	printf("jj=%d, first=%d, last=%d\n", jj, gr.tabmagalfaRange[jj].first, gr.tabmagalfaRange[jj].last);
//            }

            gr.IMMAXalfa = im;

            if(rank == 0)
                printf("gr.IMMAXalfa = %d\n", gr.IMMAXalfa);

            if(gr.IMMAXalfa >= MAGROZ -1) {
                if(rank == 0)
                    printf("gr.IMMAXalfa > MAGROZ\n");
                ProgramStop(0);
                return;
            }


            //       getc(stdin);

            for (ii=0; ii<=gr.IMMAXalfa; ii++) {

                double Omegamag = pow(Graw*konf.Mass/pow(gr.tabmagalfa[ii].rmag,3.0), 0.5);

                if(konf.flicker_dyn) {
                    gr.tabmagalfa[ii].taudynmag = 1.0/Omegamag;
                    printf("ii=%d, taudynmag=%E, oczko=%d\n", ii, gr.tabmagalfa[ii].taudynmag, gr.tabmagalfa[ii].oczko);
                }

                if(konf.flicker_visc) {
                    gr.tabmagalfa[ii].taudynmag = 1.0/Omegamag/konf.Alpha*pow(gr.tab[jr].r*Rschw/gr.tab[jr].Rmagalfa, 2.0);

                    // gr.tabmagalfa[ii].taudynmag = 1.0/Omegamag/newtm->Alpha22[jr]*pow(gr.tab[jr].r*Rschw/newtm->hd_ster[jr], 2.0);

                    printf("ii=%d, taumag_visc=%E, oczko=%d\n", ii, gr.tabmagalfa[ii].taudynmag, gr.tabmagalfa[ii].oczko);
                }
                gr.tabmagalfa[ii].lasttimedyn = 0.0;
                gr.tabmagalfa[ii].urand = konf.Alpha;
		jj = gr.tabmagalfa[ii].oczko;
		//                gr.tabmagalfa[ii].P1 = (1.0 - pow((3.0/gr.tab[jj].r/Rschw), 0.5))/(pow(Graw*konf.Mass/pow(gr.tabmagalfa[ii].rmag,3.0), 0.5));

gr.tabmagalfa[ii].P1 = (1.0 - pow((3.0*Rschw/gr.tabmagalfa[ii].rmag), 0.5))/(pow(Graw*konf.Mass/pow(gr.tabmagalfa[ii].rmag,3.0), 0.5));

 printf("ii=%d, nawias_stary=%E, nawias_nowy=%E, promien_mag[Rs]=%E\n", ii, (1.0 - pow((3.0/gr.tab[jj].r/Rschw), 0.5)), (1.0 - pow((3.0*Rschw/gr.tabmagalfa[ii].rmag), 0.5)), gr.tabmagalfa[ii].rmag/Rschw);

            }

        }

    }
    else {
        printf("Constant alfa, no flickering\n");
    }
    /////////////////////////////////////////////

    if(rank == 0)
        printf("Calculating S-curves...\n");



    if(konf.ji == 1)
    {
        scurve(&gr.tab[0]);
        WriteOneKrzywaS(0);
    }

    for (jr = konf.ji; jr <= konf.jo; jr++) {
        scurve(&gr.tab[jr]);
        WriteOneKrzywaS(jr);
    }


    scurve(&gr.tab[PROZ-1]);
    WriteOneKrzywaS(PROZ-1);

    //  scurve(&gr.tab[3]);


    for (jr = konf.ji; jr <= konf.jo; jr++)        //initial values: lower branch of the S-curve. Point 20 (mdot) arbitrarily chosen
    {
        curtm->sig_ster[jr]=gr.tab[jr].sigma[20];
        curtm->temp_ster[jr]=gr.tab[jr].temper[20];
        curtm->Teff[jr]=gr.tab[jr].Tef[20];

        curtm->hd_ster[jr]=gr.tab[jr].h[20];
        curtm->nu_ster[jr]=0.66*konf.Alpha/gr.tab[jr].kep;
        curtm->nu_ster[jr]=curtm->nu_ster[jr]
                           *(curtm->hd_ster[jr]/curtm->sig_ster[jr]*4/3*Sigmabol/Crad*pow(curtm->temp_ster[jr], 4.0) +
                             Kbol/Mprot*curtm->temp_ster[jr]);
        curtm->signu[jr]=curtm->sig_ster[jr]*curtm->nu_ster[jr];
        curtm->s[jr]=curtm->sig_ster[jr]*gr.y[jr];
        curtm->snu[jr]=curtm->s[jr]*curtm->nu_ster[jr];


        curtm->Ptot[jr]=Kbol/Mprot*curtm->sig_ster[jr]/curtm->hd_ster[jr]*curtm->temp_ster[jr]+4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr],4.0);
        ksi=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr],3.0)*Mprot/Kbol*curtm->hd_ster[jr]/curtm->sig_ster[jr];
        curtm->beta[jr]=1.0/(1.0+ksi);

        curtm->Pgas[jr]=Kbol/Mprot*curtm->sig_ster[jr]/curtm->hd_ster[jr]*curtm->temp_ster[jr];
        curtm->Prad[jr]=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr],4.0);

        //hd_old[jr]=Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/curtm->sig_ster[jr];
        curtm->hd_old[jr]=curtm->hd_ster[jr];



        if(konf.wymiana_masy) {
 //           curtm->sig_cor[jr] = 0.1/konf.Kappa; // glebokosc optyczna poczatkowa korony std 1.0
              curtm->sig_cor[jr] = 1.0/konf.Kappa;

            curtm->s_cor[jr] = curtm->sig_cor[jr]*gr.y[jr];
            curtm->h_cor[jr] = gr.tab[jr].r*Rschw;            								//Thick corona H=r
            curtm->t_cor[jr] = Graw*Mprot/Kbol*konf.Mass/curtm->h_cor[jr]; 	 //virial temperature of corona
            curtm->nu_cor[jr]=0.66*konf.AlphaCor/gr.tab[jr].kep; 							//Viscosity in corona
            curtm->nu_cor[jr]=curtm->nu_cor[jr]*Kbol/Mprot*curtm->t_cor[jr]; //P(corona)=Pgas
            curtm->snu_cor[jr] = curtm->s_cor[jr]*curtm->nu_cor[jr];
            curtm->lasttimedynprom[jr] = 0.0;
        }

    }

    //bounary conditions

    if(konf.ji == 1) {

        curtm->sig_ster[0]=0.0;                         //inner edge: T=Sigma = 0
        curtm->temp_ster[0]=0.0;
        curtm->signu[0]=0.0;
        curtm->s[0]=0.0;
        curtm->snu[0]=0.0;
        curtm->Teff[0]=0.0;
        curtm->flux[0]=0.0;
        //    curtm->vr[0]=0.0;
        curtm->curdt = konf.DtStart;


        curtm->vr_dysk[0]= 6.0/curtm->s[1]/gr.y[1]*(curtm->snu[1]-curtm->snu[0])/gr.dy[1];


        if(konf.wymiana_masy) {

            curtm->sig_cor[0]=1.0/konf.Kappa;
            //curtm->sig_cor[0]=0.1/konf.Kappa; // glebokosc opt. korony na r_in;

            curtm->t_cor[0] = Graw*Mprot/Kbol*konf.Mass/gr.tab[0].r/Rschw;
            curtm->s_cor[0] = curtm->sig_cor[0]*gr.y[0];
            curtm->nu_cor[0]=0.66*konf.AlphaCor/gr.tab[0].kep*Kbol/Mprot*curtm->t_cor[0];
            curtm->snu_cor[0] = curtm->s_cor[0]* curtm->nu_cor[0];
        }

        curtm->flux[0] = 0.0;
        curtm->flux2[0] = 0.0;
        curtm->flux3[0] = 0.0;
        curtm->fluxcor[0]=0.0;

        curtm->urand[0] = 0.0;
        curtm->urandalfa[0] = 0.0;
        curtm->Alpha22[0]=konf.Alpha;

        setRandomTable(&curtm->randd[0], 12345*(0+1),65435*(0+1),34221*(0+1),43543*(0+1));
    } //konf.ji

    //outer edge: fixed mdot


    curtm->Teff[PROZ-1]=3.0/8.0/Pi/Sigmabol*Graw*konf.Mass*Mext/pow((gr.tab[PROZ-1].r)*Rschw, 3.0)*(1.0-sqrt(3.0/gr.tab[PROZ-1].r));
    curtm->Teff[PROZ-1]=pow(curtm->Teff[PROZ-1], 0.25);

    //wyszukiwanie sigma
    for(j=0; j < MROZ; j++) {
        if(gr.tab[PROZ-1].Tef[j] >= curtm->Teff[PROZ-1])
            break;
    }
    //getc(stdin);
    if(j == MROZ) {
        //if(rank == 0)
        printf("Sigma_out not found for Mext=%E\n", Mext);
    }
    //if(rank == 0)
    printf("Sigma_out found for Mext=%E, at j=%d\n", Mext, j);


    for (jr = konf.ji; jr <= konf.jo; jr++)        //initial values: lower branch of the S-curve. Point (mdot) from outer condition instead of 20
    {
        curtm->sig_ster[jr]=gr.tab[jr].sigma[j];
        curtm->temp_ster[jr]=gr.tab[jr].temper[j];
        curtm->Teff[jr]=gr.tab[jr].Tef[j];

        curtm->hd_ster[jr]=gr.tab[jr].h[j];
        curtm->nu_ster[jr]=0.66*konf.Alpha/gr.tab[jr].kep;
        curtm->nu_ster[jr]=curtm->nu_ster[jr]
                           *(curtm->hd_ster[jr]/curtm->sig_ster[jr]*4/3*Sigmabol/Crad*pow(curtm->temp_ster[jr], 4.0) +
                             Kbol/Mprot*curtm->temp_ster[jr]);
        curtm->signu[jr]=curtm->sig_ster[jr]*curtm->nu_ster[jr];
        curtm->s[jr]=curtm->sig_ster[jr]*gr.y[jr];
        curtm->snu[jr]=curtm->s[jr]*curtm->nu_ster[jr];


        curtm->Ptot[jr]=Kbol/Mprot*curtm->sig_ster[jr]/curtm->hd_ster[jr]*curtm->temp_ster[jr]+4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr],4.0);
        ksi=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr],3.0)*Mprot/Kbol*curtm->hd_ster[jr]/curtm->sig_ster[jr];
        curtm->beta[jr]=1.0/(1.0+ksi);

        curtm->Pgas[jr]=Kbol/Mprot*curtm->sig_ster[jr]/curtm->hd_ster[jr]*curtm->temp_ster[jr];
        curtm->Prad[jr]=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[jr],4.0);

        //hd_old[jr]=Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/curtm->sig_ster[jr];
        curtm->hd_old[jr]=curtm->hd_ster[jr];
        curtm->vr_dysk[jr]= 6.0/curtm->s[jr]/gr.y[jr]*(curtm->snu[jr]-curtm->snu[jr])/gr.dy[jr];


        if(konf.wymiana_masy) {
            //curtm->sig_cor[jr] = 0.1/konf.Kappa; // glebokosc optyczna poczatkowa korony std 1.0
            curtm->sig_cor[jr] = 1.0/konf.Kappa;

            curtm->s_cor[jr] = curtm->sig_cor[jr]*gr.y[jr];
            curtm->h_cor[jr] = gr.tab[jr].r*Rschw;            								//Thick corona H=r
            curtm->t_cor[jr] = Graw*Mprot/Kbol*konf.Mass/curtm->h_cor[jr]; 	 //virial temperature of corona
            curtm->nu_cor[jr]=0.66*konf.AlphaCor/gr.tab[jr].kep; 							//Viscosity in corona
            curtm->nu_cor[jr]=curtm->nu_cor[jr]*Kbol/Mprot*curtm->t_cor[jr]; //P(corona)=Pgas
            curtm->snu_cor[jr] = curtm->s_cor[jr]*curtm->nu_cor[jr];
            curtm->lasttimedynprom[jr] = 0.0;
        }

        curtm->urand[jr] = 0.0;
        curtm->urandalfa[jr] = 0.0;
        curtm->Alpha22[jr]=konf.Alpha;

        setRandomTable(&curtm->randd[jr], 12345*(jr+1),65435*(jr+1),34221*(jr+1),43543*(jr+1));
    }

    if(konf.jo == PROZ-2) {

        curtm->sig_ster[PROZ-1]=gr.tab[PROZ-1].sigma[j];
        curtm->temp_ster[PROZ-1]=gr.tab[PROZ-1].temper[j];

        curtm->s[PROZ-1]=curtm->sig_ster[PROZ-1]*gr.y[PROZ-1];

        curtm->hd_ster[PROZ-1]=gr.tab[PROZ-1].h[j];
        curtm->nu_ster[PROZ-1]=0.66*konf.Alpha/gr.tab[PROZ-1].kep;
        curtm->nu_ster[PROZ-1]=curtm->nu_ster[PROZ-1]*(curtm->hd_ster[PROZ-1]/curtm->sig_ster[PROZ-1]*4/3*Sigmabol/Crad
                               *pow(curtm->temp_ster[PROZ-1], 4.0) +
                               Kbol/Mprot*curtm->temp_ster[PROZ-1]);
        curtm->snu[PROZ-1]=curtm->s[PROZ-1]*curtm->nu_ster[PROZ-1];

        curtm->hd_ster[PROZ-1]=gr.tab[PROZ-1].h[j];
        curtm->Ptot[PROZ-1]=Kbol/Mprot*curtm->sig_ster[PROZ-1]/curtm->hd_ster[PROZ-1]
                            *curtm->temp_ster[PROZ-1]+4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[PROZ-1],4.0);

        ksi=4.0/3.0*Sigmabol/Crad*pow(curtm->temp_ster[PROZ-1],3.0)*Mprot/Kbol*curtm->hd_ster[PROZ-1]/curtm->sig_ster[PROZ-1];
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
        curtm->urandalfa[PROZ-1] = 0.0;
        curtm->Alpha22[PROZ-1]=konf.Alpha;

        setRandomTable(&curtm->randd[PROZ-1], 12345*(PROZ),65435*(PROZ),34221*(PROZ),43543*(PROZ));

    } //konf.jo


    if(size > 1)
        bval1(curtm);                   //Wypelnienie wartosci brzegowych snu i snu_cor, do pierwszego wywolania derivatives

}




/**************************************************************/
void ProgramStop(int code)
/**************************************************************/
{
    //    signal(SIGINT, SIG_IGN);
    //    signal(SIGTERM, SIG_IGN);
    //    signal(SIGQUIT, SIG_IGN);
    //    signal(SIGHUP, SIG_IGN);

    if(rank == 0)
        printf("STOP BEGIN:\n");


    FinalizeOutput();

    if(rank == size - 1)
        DeInitOuterBoundary();


#ifdef MPI_USED
    MPI_Finalize( );
#endif

    if(rank == 0)
        printf("STOP END\n");

    exit(code);
}




/**************************************************************/
void derivatives(double Dt, SOLINTIME *curtm, DERIVAT *ddt)
//void derivatives(double Dt, DERIVAT *ddt)
//Brzeg wypelniony w Calculations po wypelnieniu bval1
/**************************************************************/
{
    register int jr;

    double spy[PROZ];
    double d_spy_dr1[PROZ];
    double d_temp_dr1[PROZ];
    double d_snu_dr2[PROZ];
    double d_snu_cor_dr2[PROZ];

    for (jr = konf.ji-1; jr <= konf.jo+1; jr++)
        spy[jr] = curtm->s[jr]/gr.y[jr];

    DerivDr1(spy, curtm->vr_dysk, d_spy_dr1);
    DerivDr1(curtm->temp_ster, curtm->vr_dysk, d_temp_dr1);

    DerivDr2(curtm->snu, d_snu_dr2);

    if(konf.wymiana_masy && (!konf.PBK || (konf.PBK && (ptime > konf.PBK_time)))) 
        DerivDr2(curtm->snu_cor, d_snu_cor_dr2);


    for (jr = konf.ji; jr <= konf.jo; jr++) {

        /*Diffusion equation for Sigma x nu*/

//        ddt->dsdt[jr]=12.0/gr.y[jr]/gr.y[jr]*(curtm->snu[jr+1]+curtm->snu[jr-1]-2.0*curtm->snu[jr])/gr.dy[jr]/gr.dy[jr] - curtm->mz[jr]*gr.y[jr];
        ddt->dsdt[jr]=12.0/gr.y[jr]/gr.y[jr]*d_snu_dr2[jr];

		if(konf.wymiana_masy && (!konf.PBK || (konf.PBK && (ptime > konf.PBK_time))))
			ddt->dsdt[jr] -= curtm->mz[jr]*gr.y[jr];

        //	printf("DER jr=%d,y=%E  %E %E %E %E\n", jr, gr.y[jr], ddt->dsdt[jr], d_snu_dr2[jr], curtm->mz[jr], curtm->snu[jr]);

        if(konf.wymiana_masy) {

            if (!konf.PBK || (konf.PBK && (ptime > konf.PBK_time))) {

//                ddt->dsdt_cor[jr]=12.0/gr.y[jr]/gr.y[jr]*(curtm->snu_cor[jr+1]+curtm->snu_cor[jr-1]-2.0*curtm->snu_cor[jr])/gr.dy[jr]/gr.dy[jr]
//                                  + curtm->mz[jr]*gr.y[jr];
                ddt->dsdt_cor[jr]=12.0/gr.y[jr]/gr.y[jr]*d_snu_cor_dr2[jr] + curtm->mz[jr]*gr.y[jr];

            }
            else {

                ddt->dsdt_cor[jr]= 0.0;

            }

        }


        /*Energy equation for Tc*/

        ddt->dtempdt[jr] = (4.0-3.0*curtm->beta[jr])/(12.0-10.5*curtm->beta[jr])*curtm->temp_ster[jr]*
                           (ddt->dsdt[jr]/gr.y[jr]/curtm->sig_ster[jr]
                            //- (curtm->hd_ster[jr]-curtm->hd_old[jr])/*/curtm->curdt*//curtm->hd_ster[jr]
                            //- 12.0/curtm->s[jr]/y[jr]/y[jr]/curtm->hd_ster[jr]*(curtm->snu[jr+1]-curtm->snu[jr-1])/2.0/dy[jr]*
                            //(curtm->hd_ster[jr+1]-curtm->hd_ster[jr-1])/2.0/y[jr]
                            //+ 12.0/curtm->s[jr]/curtm->s[jr]*(curtm->s[jr+1]/y[jr+1]-curtm->s[jr-1]/y[jr-1])/2.0/dy[jr]/y[jr]*(curtm->snu[jr+1]-curtm->snu[jr-1])/2.0/dy[jr]);
//                            + curtm->vr_dysk[jr]*2.0/gr.y[jr]/curtm->sig_ster[jr]*(curtm->s[jr+1]/gr.y[jr+1]-curtm->s[jr-1]/gr.y[jr-1])/2.0/gr.dy[jr]);
                            + curtm->vr_dysk[jr]*2.0/gr.y[jr]/curtm->sig_ster[jr]*d_spy_dr1[jr]);


        /*
        if(heating)
        {
            ddt->dtempdt[jr] = ddt->dtempdt[jr] +
        ((1.0-Xicor)*(1.0-etajet3)*C1*1.5*Alpha22[jr]*gr.tab[jr].kep*curtm->temp_ster[jr]/(12.0-10.5*curtm->beta[jr])/C6 -
        C2*4.0/(12.0-10.5*curtm->beta[jr])*Sigmabol/3.0/Kappa*pow(curtm->temp_ster[jr],5.0)
        /curtm->Ptot[jr]/curtm->hd_ster[jr]*y[jr]/(curtm->s[jr])/C6);

        }
        else
        {

        */

        /*
               ddt->dtempdt[jr] = ddt->dtempdt[jr] + (C1*1.5*curtm->Alpha22[jr]*gr.tab[jr].kep*curtm->temp_ster[jr]/(12.0-10.5*curtm->beta[jr])/C6 -
                                                      C2/((1.0-Xicor)*(1.0-etajet3))*4.0/(12.0-10.5*curtm->beta[jr])*Sigmabol/3.0/Kappa*pow(curtm->temp_ster[jr],5.0)
                                                      *y[jr]/(curtm->s[jr])/(curtm->Ptot[jr]*curtm->hd_ster[jr])/C6);
        				       */


        ddt->dtempdt[jr] += (curtm->Qplus[jr]-curtm->Qminus[jr])*curtm->temp_ster[jr]/(12.0-10.5*curtm->beta[jr])/konf.C6/(curtm->Ptot[jr]*curtm->hd_ster[jr]);

//        ddt->dtempdt[jr] -= curtm->vr_dysk[jr]*2.0/gr.y[jr]*(curtm->temp_ster[jr+1]-curtm->temp_ster[jr-1])/2.0/gr.dy[jr];
        ddt->dtempdt[jr] -= curtm->vr_dysk[jr]*2.0/gr.y[jr]*d_temp_dr1[jr];

//ddt->dtempdt[jr] =0.0;
//ddt->dsdt[jr] = 0.0;
//ddt->dsdt_cor[jr]= 0.0;
//ddt->dsdt[jr]=curtm->s[jr]/100000000000.0;
    }

    //getc(stdin);

}

/**************************************************************/
//void calculations(double Dt, DERIVAT *ddt, double derivscale, int calclum)
void calculations(double Dt, SOLINTIME *curtm, SOLINTIME *newtm, DERIVAT *ddt, double derivscale, int calclum)

/**************************************************************/
{
    register int jr;
    double dtemp;

    //Na paczatku przyjmujemy, ze nowe bedzie sie rownalo biezacemu
    memcpy(newtm, curtm, sizeof(SOLINTIME));

    for (jr = konf.ji; jr <= konf.jo; jr++)   //1. petla
    {

        newtm->s[jr] = curtm->s[jr]+ddt->dsdt[jr]*derivscale;

        if(newtm->s[jr] < 1.e-20)
           newtm->s[jr] = 1.e-20;


        newtm->sig_ster[jr]=newtm->s[jr]/gr.y[jr];

        //	printf("jr=%d, %E %E %E %E %E\n", jr,  newtm->sig_ster[jr], newtm->s[jr], gr.y[jr], curtm->s[jr], ddt->dsdt[jr]*derivscale);

        newtm->temp_ster[jr]=curtm->temp_ster[jr]+ddt->dtempdt[jr]*derivscale;

        if(newtm->temp_ster[jr] < 1.e-7)
           newtm->temp_ster[jr] = 1.e-7;

        dtemp = Kbol/Mprot*newtm->sig_ster[jr]/curtm->hd_ster[jr]*newtm->temp_ster[jr];

        newtm->Ptot[jr]= dtemp + 4.0/3.0*Sigmabol/Crad*pow(newtm->temp_ster[jr], 4.0);
        newtm->Pgas[jr]= dtemp;
        newtm->Prad[jr]=4.0/3.0*Sigmabol/Crad*pow(newtm->temp_ster[jr], 4.0);


        newtm->hd_old[jr]=curtm->hd_ster[jr];
        newtm->hd_ster[jr]=newtm->Ptot[jr]/gr.tab[jr].kep/gr.tab[jr].kep/newtm->sig_ster[jr]/sqrt(konf.C3);

        if(konf.ipsqrt) {
            newtm->nu_ster[jr]=sqrt(newtm->Pgas[jr]*newtm->Ptot[jr])*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*newtm->hd_ster[jr]/newtm->sig_ster[jr];

        } else {
            if (konf.iptot) {
                newtm->nu_ster[jr]=newtm->Ptot[jr]*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*newtm->hd_ster[jr]/newtm->sig_ster[jr];
            } else {
                newtm->nu_ster[jr]=newtm->Pgas[jr]*0.66*curtm->Alpha22[jr]/gr.tab[jr].kep*newtm->hd_ster[jr]/newtm->sig_ster[jr];
            }

        }

        if(konf.Dzet != 0.0) {

            //	double etajet;
            etajet = 1.0 - 1.0/(1.0+ konf.Dzet*pow(2*Pi*gr.tab[jr].r*Rschw*newtm->sig_ster[jr]*curtm->vr_dysk[jr]*Crad*Crad/16.0/Ledd  , 2.0));
            //printf("jr=%d, etajet = %E %E %E\n", jr, etajet, newtm->sig_ster[jr], curtm->vr_dysk[jr]);
        }
        else {
            etajet = konf.etajet3;
        }

        newtm->Qminus[jr]=Crad*konf.C2/(1.0-konf.Xicor)/(1.0-etajet)*newtm->Prad[jr]/konf.Kappa/newtm->sig_ster[jr];

        //        newtm->Qminus[jr]=Crad*konf.C2/(1.0-konf.Xicor)/(1.0-konf.etajet3)*newtm->Prad[jr]/konf.Kappa/newtm->sig_ster[jr];


        newtm->snu[jr]=newtm->s[jr]*newtm->nu_ster[jr];


        double ksi;
        ksi=4.0/3.0*Sigmabol/Crad*pow(newtm->temp_ster[jr], 3.0)*Mprot/Kbol*newtm->hd_ster[jr]/newtm->sig_ster[jr];
        newtm->beta[jr]=1.0/(1.0+ksi);

	    if(konf.alfa_flickering)
    	    newtm->urandalfa[jr] = CalcUrandAlfa(jr, curtm->urandalfa[jr], ptime+Dt, &newtm->randd[jr], newtm);
    	else    
    	    newtm->urandalfa[jr] = konf.Alpha;

        if(konf.viscosity) {
	  newtm->Alpha22[jr]=Flicker(jr, newtm->urandalfa[jr], gr.tab[jr].r*Rschw, ptime+Dt, curtm->hd_ster[jr])*(1.0+ksi/konf.xi0)/(1.0+pow((ksi/konf.xi0), 2.0));
        } else {
	  newtm->Alpha22[jr]=Flicker(jr, newtm->urandalfa[jr], gr.tab[jr].r*Rschw, ptime+Dt, curtm->hd_ster[jr]);
        }

        if(konf.ipsqrt) {
            newtm->Qplus[jr]=konf.C1*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->hd_ster[jr]*sqrt(newtm->Pgas[jr]*newtm->Ptot[jr]);

            //            newtm->Qplus[jr]=konf.C1*(1.0-konf.Xicor)*(1.0-konf.etajet3)*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->hd_ster[jr]*sqrt(newtm->Pgas[jr]*newtm->Ptot[jr]);
        } else {

            if (konf.iptot) {
                newtm->Qplus[jr]=konf.C1*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->hd_ster[jr]*newtm->Ptot[jr];
                //                newtm->Qplus[jr]=konf.C1*(1.0-konf.Xicor)*(1.0-konf.etajet3)*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->hd_ster[jr]*newtm->Ptot[jr];
            } else {
                newtm->Qplus[jr]=konf.C1*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->hd_ster[jr]*newtm->Pgas[jr];
                //                newtm->Qplus[jr]=konf.C1*(1.0-konf.Xicor)*(1.0-konf.etajet3)*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->hd_ster[jr]*newtm->Pgas[jr];
            }

        }

        newtm->Qminus[jr]=Crad*konf.C2/(1.0-konf.Xicor)/(1.0-etajet)*newtm->Prad[jr]/konf.Kappa/newtm->sig_ster[jr];
        //        newtm->Qminus[jr]=Crad*konf.C2/(1.0-konf.Xicor)/(1.0-konf.etajet3)*newtm->Prad[jr]/konf.Kappa/newtm->sig_ster[jr];


        if(konf.wymiana_masy) {


            newtm->s_cor[jr] = curtm->s_cor[jr]+ddt->dsdt_cor[jr]*derivscale;

            newtm->sig_cor[jr] = newtm->s_cor[jr]/gr.y[jr];

            newtm->nu_cor[jr]=0.66*konf.AlphaCor/gr.tab[jr].kep;

            newtm->nu_cor[jr]=newtm->nu_cor[jr]*Kbol/Mprot*curtm->t_cor[jr];

            newtm->snu_cor[jr] = newtm->s_cor[jr]*newtm->nu_cor[jr];

        }
    } //1. petla



    if(size > 1)
        bval1(newtm);                   //Wypelnienie wartosci brzegowych snu i snu_cor

//    printf("N %d rank: %d %2d %2d %E %E\n", ntime, rank, konf.ji, konf.jo, newtm->s[konf.ji], newtm->s[konf.jo]);
//    printf("B %d rank: %d %2d %2d %E %E\n", ntime, rank, konf.ji-1, konf.jo+1, newtm->s[konf.ji-1], newtm->s[konf.jo+1]);


    double d_snu_dr1[PROZ];
    double d_snu_cor_dr1[PROZ];

    DerivDr1(newtm->snu, curtm->vr_dysk, d_snu_dr1);

    if(konf.wymiana_masy)
        DerivDr1(newtm->snu_cor, curtm->vr_cor, d_snu_cor_dr1);


    for (jr = konf.ji; jr <= konf.jo; jr++)   //1.1 petla
    {
//        newtm->vr_dysk[jr]= 6.0/newtm->s[jr]/gr.y[jr]*(newtm->snu[jr+1]-newtm->snu[jr-1])/2.0/gr.dy[jr];
        newtm->vr_dysk[jr]= 6.0/newtm->s[jr]/gr.y[jr]*d_snu_dr1[jr];

        if(konf.wymiana_masy) {

//            newtm->vr_cor[jr] =  6.0/newtm->s_cor[jr]/gr.y[jr]*(newtm->snu_cor[jr+1]-newtm->snu_cor[jr-1])/2.0/gr.dy[jr];
            newtm->vr_cor[jr] =  6.0/newtm->s_cor[jr]/gr.y[jr]*d_snu_cor_dr1[jr];

            //Mass exchange: various prescriptions:

            if (!konf.PBK || (konf.PBK && (ptime > konf.PBK_time))) {

                if(konf.wymiana_flux) {

		  // prescription copied directly from old 2005 paper code version, according to Eq. (17):

		  /* newtm->mz[jr] = 3.0/4.0/(gr.tab[jr].r*Rschw)*
                    (1.-pow((3./gr.tab[jr].r),1.5)*(gr.tab[jr].r-1.)/2.)*
                    (0.5*konf.B1*newtm->sig_cor[jr]*newtm->vr_cor[jr]+konf.B2*newtm->sig_ster[jr]*newtm->vr_dysk[jr]);
		  */

		  //prescription according to Eq. (26)
		  
		 newtm->mz[jr] = 3.0/2.0*gr.tab[jr].kep*(konf.B2*newtm->Alpha22[jr]*newtm->Ptot[jr]*newtm->hd_ster[jr]+konf.B1*0.5*konf.AlphaCor*newtm->sig_cor[jr]*Graw*konf.Mass/gr.tab[jr].r/Rschw)/Graw/konf.Mass*gr.tab[jr].r*Rschw*(1-pow((3./gr.tab[jr].r),1.5)*(gr.tab[jr].r-1.)/2.);

/* ====================================================== */

		  //         newtm->mz[jr] = 3.0/2.0*gr.tab[jr].kep*(konf.B2*newtm->Alpha22[jr]*newtm->Ptot[jr]*newtm->hd_ster[jr]+
                  //                                          konf.B1*0.5*konf.AlphaCor*newtm->sig_cor[jr]*Graw*konf.Mass/gr.tab[jr].r/Rschw)/Graw/konf.Mass*gr.tab[jr].r*Rschw* /*3.0/4.0/(gr.tab[jr].r*Rschw)* */
		  // (1-pow((3./gr.tab[jr].r),1.5)*(gr.tab[jr].r-1.)/2.); //if bracket with 3/4/r/Rs commented out -> timestep goes extremely down; if not -> very small corona luminosiy (check that also with Sigma_cor = 1.0/Kappa)

                    /*3.0/4.0/(gr.tab[jr].r*Rschw)*
                    (1-pow((3./gr.tab[jr].r),1.5)*(gr.tab[jr].r-1.)/2.)*
                    (0.5*konf.B1*newtm->sig_cor[jr]*newtm->vr_cor[jr]+konf.B2*newtm->sig_ster[jr]*newtm->vr_dysk[jr]);
                    */

                } else {

                    newtm->BzMax[jr] = konf.BetaS*sqrt(4.0*Pi*newtm->Alpha22[jr]*newtm->Ptot[jr]);

                    if(konf.grid_radius) {

                        if (ptime + Dt >= curtm->lasttimedynprom[jr] + gr.tab[jr].taudyn) {
                            //    printf("jest markoff jr=%d time = %E last=%E tau=%E\n", jr, ptime+Dt, curtm->lasttimedynprom[jr], gr.tab[jr].taudyn);

                            newtm->urand[jr] = Markoff(curtm->urand[jr], &newtm->randd[jr]);
                            newtm->lasttimedynprom[jr] = ptime + Dt;

                            //   printf("cur=%E new=%E\n", curtm->urand[jr], newtm->urand[jr]);
                            //   getc(stdin);

                        }
                        else {

                            //  printf("nie jest markoff jr=%d time = %E last=%E tau =%E\n", jr, ptime+Dt, curtm->lasttimedynprom[jr], gr.tab[jr].taudyn);

                            newtm->urand[jr] = curtm->urand[jr];
                        }


                    } else {

                        int jm, ile,ilemarkof ;
                        double usum, um;

                        usum = 0.0;
                        ile = 0;
                        ilemarkof=0;

                        for (jm=0; jm<=gr.IMMAX; jm++) {

                            if(gr.tabmag[jm].oczko == jr) {

			      // gr.tabmagalfa[ii].taudynmag = 1.0/Omegamag/newtm->Alpha22[jr]*pow(gr.tab[jr].r*Rschw/newtm->hd_ster[jr], 2.0);

                                if (ptime + Dt >= gr.tabmag[jm].lasttimedyn + gr.tabmag[jm].taudynmag) {

                                    //	  printf("jest markoff jr=%d jm=%d time = %E last=%E taudyn=%E taumag=%E\n", jr, jm, ptime+Dt, gr.tabmag[jm].lasttimedyn, gr.tab[jr].taudyn, gr.tabmag[jm].taudynmag);
                                    //	getc(stdin);

                                    gr.tabmag[jm].lasttimedyn = ptime + Dt;

                                    um = Markoff(gr.tabmag[jm].urand, &newtm->randd[jr]);
                                    gr.tabmag[jm].urand = um;

                                    ilemarkof++;
                                }
                                else {
                                    um = gr.tabmag[jm].urand;

                                }
                                usum = usum + um;
                                ile++;

                            } else {
                                //	newtm->urand[jr] = curtm->urand[jr];
                                continue;

                            }

                        }


                        if(ile>0)
                        {
                            newtm->urand[jr] = usum/ile;

                        } else
                            newtm->urand[jr] = curtm->urand[jr];

                    } //else


                    newtm->Bz[jr] = newtm->BzMax[jr]*newtm->urand[jr];

                    //      newtm->mz[jr] = newtm->Bz[jr]*newtm->Bz[jr]/8.0/Pi/Crad;
                    //		if(newtm->mz[jr]!=0.0)		printf("jr = %d, Old mz = %E\n", jr, newtm->mz[jr]);

                    newtm->Valvf[jr] = sqrt(newtm->Alpha22[jr])*gr.tab[jr].kep*newtm->hd_ster[jr];

                    switch(konf.case_mz) {

                    case 'A':
                        newtm->mz[jr] = newtm->Bz[jr]*newtm->Bz[jr]/8.0/Pi*newtm->Valvf[jr]/(Graw*konf.Mass/gr.tab[jr].r/Rschw);
                        break;

                    case 'B':
                        newtm->mz[jr] = newtm->Bz[jr]*newtm->Bz[jr]/8.0/Pi*Crad/(Graw*konf.Mass/gr.tab[jr].r/Rschw);
                        break;
                    case 'C':
                        newtm->mz[jr] = newtm->Bz[jr]*newtm->Bz[jr]/8.0/Pi*newtm->Valvf[jr]/(Graw*konf.Mass/gr.tab[jr].r/Rschw);
                        newtm->mz[jr] /= (1.0-sqrt(3.0/gr.tab[jr].r));
                        //(1.0-pow((3.0/gr.tab[jr].r),1.5)*(gr.tab[jr].r-1.)/2.);
                        break;
                    default:
                        if(rank == 0)
                            printf("Case not defined case_mz = '%c'\n", konf.case_mz);
                        break;
                    }

                }

            }
            else {
                newtm->BzMax[jr] = 0.0;
                newtm->urand[jr] = curtm->urand[jr];

                newtm->Bz[jr] = 0.0;
                newtm->Valvf[jr] =0.0;
                newtm->mz[jr] = 0.0;
            }
        }
    }    //2. petla


    if(calclum) {

        newtm->lum=0.0;
        newtm->lum2=0.0;
        newtm->lum3=0.0;
        newtm->lumcor=0.0;

        double scorpy[PROZ];
        double d_scorpy_dr1[PROZ];

        for (jr = konf.ji-1; jr <= konf.jo+1; jr++)
            scorpy[jr] = newtm->s_cor[jr]/gr.y[jr];

        DerivDr1(scorpy, newtm->vr_cor, d_scorpy_dr1);


        for (jr = konf.ji; jr <= konf.jo; jr++) //2. petla
        {

            newtm->Teff[jr]=newtm->temp_ster[jr]/pow(0.75*(newtm->sig_ster[jr]*konf.Kappa+0.66), 0.25)*konf.C5;

            newtm->logflux=4*log10(newtm->Teff[jr])+log10(Sigmabol);

            /* Local flux*/

            newtm->flux[jr]=pow(10, newtm->logflux);


            //            newtm->flux2[jr] =  konf.C2/((1.0-konf.Xicor)*(1.0-konf.etajet3))*4.0*Sigmabol/3.0/konf.Kappa*pow(newtm->temp_ster[jr],4.0)
            //                               *gr.y[jr]/(newtm->s[jr]);

            newtm->flux2[jr] =  konf.C2*4.0*Sigmabol/3.0/konf.Kappa*pow(newtm->temp_ster[jr],4.0)*gr.y[jr]/(newtm->s[jr]);

            /*Frad*/

            newtm->Teff2[jr] = pow(newtm->flux2[jr]/Sigmabol, 0.25);

            newtm->flux3[jr]= konf.C1*1.5*newtm->Alpha22[jr]*gr.tab[jr].kep*newtm->Ptot[jr]*newtm->hd_ster[jr];          /*Fvisc*/

            newtm->Teff3[jr]= pow(newtm->flux3[jr]/Sigmabol, 0.25);

            if(konf.wymiana_masy) {
                newtm->fluxcor[jr]= 1.5*konf.AlphaCor*gr.tab[jr].kep*newtm->sig_cor[jr]*Graw*konf.Mass/gr.tab[jr].r/Rschw
                                    + newtm->sig_cor[jr]*Graw*konf.Mass/gr.tab[jr].r/Rschw*newtm->vr_cor[jr]*(1.5/gr.tab[jr].r/Rschw
//                                            + 1.0/newtm->s_cor[jr]*(newtm->s_cor[jr+1]/gr.y[jr+1]-newtm->s_cor[jr-1]/gr.y[jr-1])/gr.dy[jr]);
                                            + 1.0/newtm->s_cor[jr]*d_scorpy_dr1[jr]);


                if(newtm->fluxcor[jr] < 0.0)
                    newtm->fluxcor[jr]=0.0;

            }


        }

        //            }

        //trzeba wypelnic lewy brzeg: flux, flux2, flux3, fluxcor


        if(size > 1)
            bval2(newtm);                   //Wypelnienie wartosci brzegowych


        /*Fcor*/

        for (jr = konf.ji; jr <= konf.jo; jr++) //3. petla
        {

            newtm->lum += Pi*(newtm->flux[jr-1]*gr.tab[jr-1].r+newtm->flux[jr]*gr.tab[jr].r)*(gr.tab[jr].r-gr.tab[jr-1].r)*Rschw*Rschw;

            newtm->lum2 += Pi*(newtm->flux2[jr-1]*gr.tab[jr-1].r+newtm->flux2[jr]*gr.tab[jr].r)*(gr.tab[jr].r-gr.tab[jr-1].r)*Rschw*Rschw;

            newtm->lum3 += Pi*(newtm->flux3[jr-1]*gr.tab[jr-1].r+newtm->flux3[jr]*gr.tab[jr].r)*(gr.tab[jr].r-gr.tab[jr-1].r)*Rschw*Rschw;

            if(konf.wymiana_masy) {

                newtm->lumcor += Pi*(newtm->fluxcor[jr-1]*gr.tab[jr-1].r+newtm->fluxcor[jr]*gr.tab[jr].r)*(gr.tab[jr].r-gr.tab[jr-1].r)*Rschw*Rschw;
            }

        }

    }



}



/*********************************************/
double Markoff(double urand, RANDOMDEF *rd)
/*********************************************/
{

    //    double a = 0.5;
    //    double ret = urand*(-a) + (0.5-(double)rand()/RAND_MAX)*2.0;
    double ret = urand*(-0.5) + UniformRandom(rd);

    return ret;
}



/*********************************************/
static void  bval1(SOLINTIME *curtm)
//Wypelnienie wartosci brzegowych
/*********************************************/
{
#ifdef MPI_USED

    double dput1[5];
    double dget1[5];

    double dput2[5];
    double dget2[5];

    MPI_Status status[4];
    MPI_Request request[4];
    int rnum = 0;
    ///////////////////////////////////////////////////////////
    //prawy brzeg: snu, snu_cor, temp_ster, s, s_cor

    /* Send up unless I'm at the top, then receive from below */

    if (rank < size - 1) {
        dput1[0] = curtm->snu[konf.jo];
        dput1[1] = curtm->snu_cor[konf.jo];
        dput1[2] = curtm->temp_ster[konf.jo];
        dput1[3] = curtm->s[konf.jo];
        dput1[4] = curtm->s_cor[konf.jo];

        MPI_Isend( dput1, 5, MPI_DOUBLE, rank + 1, 100, MPI_COMM_WORLD, &request[rnum++]);
    }

    if (rank > 0) {
        MPI_Irecv( dget1, 5, MPI_DOUBLE, rank - 1, 100,  MPI_COMM_WORLD, &request[rnum++]);
    }

    ///////////////////////////////////////////////////////////
    //lewy brzeg: snu, snu_cor, temp_ster, s

    /* Send down unless I'm at the bottom */
    if (rank > 0) {
        dput2[0] = curtm->snu[konf.ji];
        dput2[1] = curtm->snu_cor[konf.ji];
        dput2[2] = curtm->temp_ster[konf.ji];
        dput2[3] = curtm->s[konf.ji];
        dput2[4] = curtm->s_cor[konf.ji];

        MPI_Isend( dput2, 5, MPI_DOUBLE, rank - 1, 101, MPI_COMM_WORLD, &request[rnum++]);
    }

    if (rank < size - 1) {
        MPI_Irecv( dget2, 5, MPI_DOUBLE, rank + 1, 101,  MPI_COMM_WORLD, &request[rnum++]);
    }

    MPI_Waitall(rnum, request, status);

    if (rank > 0) {
        curtm->snu[konf.ji-1]      = dget1[0];
        curtm->snu_cor[konf.ji-1]  = dget1[1];
        curtm->temp_ster[konf.ji-1]= dget1[2];
        curtm->s[konf.ji-1]        = dget1[3];
        curtm->s_cor[konf.ji-1]    = dget1[4];
    }

    if (rank < size - 1) {
        curtm->snu[konf.jo+1]      = dget2[0];
        curtm->snu_cor[konf.jo+1]  = dget2[1];
        curtm->temp_ster[konf.jo+1]= dget2[2];
        curtm->s[konf.jo+1]        = dget2[3];
        curtm->s_cor[konf.jo+1]    = dget2[4];
    }



#endif
}



/*********************************************/
static void  bval2(SOLINTIME *curtm)
//Wypelnienie wartosci brzegowych
/*********************************************/
{
#ifdef MPI_USED

    double dput[4];
    double dget[4];

    MPI_Status status[2];
    MPI_Request request[2];
    int rnum = 0;

    ///////////////////////////////////////////////////////////
    //lewy brzeg: flux, flux2, flux3, fluxcor

    /* Send up unless I'm at the top, then receive from below */
    if (rank < size - 1) {
        dput[0] = curtm->flux[konf.jo];
        dput[1] = curtm->flux2[konf.jo];
        dput[2] = curtm->flux3[konf.jo];
        dput[3] = curtm->fluxcor[konf.jo];

        MPI_Isend( dput, 4, MPI_DOUBLE, rank + 1, 104, MPI_COMM_WORLD, &request[rnum++]);
    }


    if (rank > 0) {
        MPI_Irecv( dget, 4, MPI_DOUBLE, rank - 1, 104,  MPI_COMM_WORLD, &request[rnum++]);
    }

    MPI_Waitall(rnum, request, status);

    if (rank > 0) {
        curtm->flux[konf.ji-1]  		= dget[0];
        curtm->flux2[konf.ji-1]			= dget[1];
        curtm->flux3[konf.ji-1]     = dget[2];
        curtm->fluxcor[konf.ji-1]   = dget[3];

    }

#endif
}




