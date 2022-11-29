#ifdef MPI_USED
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#ifdef HDF5_USED
#include "hdf5.h"
#endif

#include "konfig.h"
#include "dysk_zwm.h"
#include "inout.h"
#include "nsolve.h"


static FILE *sigtefiles[PROZ];
static FILE *lumfile;

static struct {
    long sigtepos[PROZ];
    long lumpos;

    double lastlum;
    double lastlumcor;

}
filepos;

void SaveFilePositions(void);
#ifdef HDF5_USED
void write_xdmf_xml(double time, int hdfNo, const char *hdf_file);
int WriteCoordsHDF(const char *FileName);
#endif

/**************************************************************/
int InitializeOutput(int restart)
/**************************************************************/
{
    char t_file[256];
    int jr;
    int minjr = (konf.ji == 1 ? 0 : konf.ji);
    int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

    const char * fattr;

    if(restart > 0)
        fattr = "r+";
    else
        fattr = "w+";


    if(konf.lsig)
        for (jr = minjr; jr <= maxjr; jr++) {
            sprintf(t_file, "sigte.%03d", jr); //files for time-dependent solutions
            if((sigtefiles[jr] = fopen(t_file, fattr))==NULL) {
                printf("Unable to open file %s\n", t_file);
                return 0;
            }
            /*
                    sprintf(t_file, "mdott.%03d", jr); //files for time-dependent solutions
                    if((mdoty[jr] = fopen(t_file,"wa"))==NULL)
                    {
                        printf("Unable to open file %s\n", t_file);
                        return 0;
                    }
            */
            if(restart > 0) {
                fseek(sigtefiles[jr], filepos.sigtepos[jr], SEEK_SET);
                ftruncate(fileno(sigtefiles[jr]), filepos.sigtepos[jr]);
            } else {
                filepos.sigtepos[jr] = 0L;
            }
        }



    if(konf.llum && rank == 0) {
        if((lumfile = fopen("lum.dat", fattr))==NULL) {
            printf("Unable to open file lum.dat\n");
            return 0;
        }

        if(restart > 0) {
            fseek(lumfile, filepos.lumpos, SEEK_SET);
            ftruncate(fileno(lumfile), filepos.lumpos);
        } else {
            filepos.lumpos = 0L;
            filepos.lastlum = 0.0;
            filepos.lastlumcor = 0.0;
        }

    }

    return 1;

}

/**************************************************************/
void SaveFilePositions(void)
/**************************************************************/
{

    if(konf.llum && rank == 0)
        filepos.lumpos = ftell(lumfile);

    if(konf.lsig) {
        int jr;
        int minjr = (konf.ji == 1 ? 0 : konf.ji);
        int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

        for (jr = minjr; jr <= maxjr; jr++) {
            filepos.sigtepos[jr] = ftell(sigtefiles[jr]);
        }
    }


}

/**************************************************************/
int WriteOutputLum(double time, int ntime, double lum2, double lum3, double lumcor)
/**************************************************************/
{
    if(konf.wymiana_masy) {
//        fprintf(lumfile, "%13lf %lf %lf %13.10lf\n",
//                time, lum2, lum3, lumcor);

      fprintf(lumfile, "%13lf %lf %13.10lf\n",
                time, lum2, lumcor);

    } else {
//        fprintf(lumfile, "%13lf %lf %lf\n",
//                time, lum2, lum3);
        fprintf(lumfile, "%13lf %lf\n",
                time, lum2);
    }
    fflush(lumfile);

    return 1;

}

/**************************************************************/
int WriteOutputSigte(double time, int ntime, SOLINTIME *newtm)
/**************************************************************/
{
    int jr;
    int minjr = (konf.ji == 1 ? 0 : konf.ji);
    int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

    for (jr = minjr; jr <= maxjr; jr++) {
        /*
                fprintf(mdoty[jr], "%E %lf\n", time, mdot_time[jr]);
                fflush(mdoty[jr]);
        */
        if(konf.wymiana_masy) {

            fprintf(sigtefiles[jr], "%E %E %E %E\n", newtm->sig_ster[jr],
                    newtm->Teff2[jr], /*newtm->Teff3[jr],*/ newtm->sig_cor[jr], newtm->mz[jr]);
        } else {

	  //if(jr == 4 || jr == 20 || jr == 45 || jr == 110 || jr == 400 ) 
//if(jr == 1 || jr == 7 || jr == 25 || jr == 50 || jr == 97 ) 
{
            fprintf(sigtefiles[jr], "%E %E %E %E %E %E \n", newtm->sig_ster[jr],
                    //newtm->Teff[jr], newtm->temp_ster[jr],
                    newtm->Teff2[jr], newtm->Ptot[jr]/*, adv_fac[jr]*/,newtm->beta[jr], /*2.0*Pi*gr.tab[jr].r*Rschw*newtm->sig_ster[jr]*newtm->vr_dysk[jr],*/ newtm->Alpha22[jr], /*1.0/newtm->Alpha22[jr]/pow(Graw*konf.Mass/pow(gr.tab[jr].r*Rschw,3.0), 0.5)*pow(gr.tab[jr].r*Rschw/newtm->hd_ster[jr], 2.0)*/ newtm->hd_ster[jr]);

 }

        }

        fflush(sigtefiles[jr]);

    }


    return 1;

}


/**************************************************************/
int WriteOutputRad(double time, int ntime, SOLINTIME *newtm)
/**************************************************************/
{

    char r5_file[64];
    FILE *plik5 = NULL;

    int jr;

    int minjr = (konf.ji == 1 ? 0 : konf.ji);
    int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

    if(rank == 0) {
        sprintf(r5_file, "rsig.%E", time);

        plik5 = fopen(r5_file,"w");
        printf("Prining radial profiles for time=%E, %d\n", time, ntime);
    }


#ifdef MPI_USED


#define ILED 11

    typedef double OUTTABLE[ILED][PROZ];

    OUTTABLE dinn;
    OUTTABLE dout;

    OUTTABLE *pdout;

    int i;

    for(i = 0; i < ILED; ++i)
        for(jr = 0; jr < minjr; ++jr)
            dinn[i][jr] = 0.0;


    for (jr = minjr; jr <= maxjr; jr++) {
        dinn[0][jr] = gr.tab[jr].r;
        dinn[1][jr] = newtm->sig_ster[jr];
        dinn[2][jr] = newtm->sig_cor[jr];
        dinn[3][jr] = newtm->mz[jr];
        dinn[4][jr] = newtm->hd_ster[jr];
        dinn[5][jr] = newtm->Ptot[jr];
        dinn[6][jr] = newtm->Pgas[jr];
        dinn[7][jr] = newtm->nu_ster[jr];
        dinn[8][jr] = newtm->vr_dysk[jr];
        dinn[9][jr] = newtm->vr_cor[jr];
        dinn[10][jr] = newtm->Alpha22[jr];
    }

    for(i = 0; i < ILED; ++i)
        for(jr = maxjr + 1; jr < PROZ; ++jr)
            dinn[i][jr] = 0.0;


    //Pobieramy sume jasnosci
    if (size > 1) {
        for(i = 0; i < ILED; ++i)
            MPI_Reduce(dinn[i], dout[i], PROZ, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        pdout = &dout;
    } else {
        pdout = &dinn;
    }

    if(rank == 0) {
        for (jr = 0; jr < PROZ; jr++)
            if(konf.wymiana_masy) {
                fprintf(plik5, "%E %E %E %E %E %E %E %E %E %E %E\n",
                        (*pdout)[0][jr],
                        (*pdout)[1][jr],
                        (*pdout)[2][jr],
                        (*pdout)[3][jr],
                        (*pdout)[4][jr],
                        (*pdout)[5][jr],
                        (*pdout)[6][jr],
                        (*pdout)[7][jr],
                        (*pdout)[8][jr],

                        (*pdout)[9][jr],

                        (*pdout)[10][jr]
                       );
            } else {

                fprintf(plik5, "%E %E %E %E %E %E %E %E\n",
                        (*pdout)[0][jr],
                        (*pdout)[1][jr],
                        (*pdout)[4][jr],

                        (*pdout)[5][jr],

                        (*pdout)[6][jr],

                        (*pdout)[7][jr],

                        (*pdout)[8][jr],

                        (*pdout)[10][jr]
                       );
            }
    }


#else

    for (jr = minjr; jr <= maxjr; jr++)
    {
        if(konf.wymiana_masy) {
            fprintf(plik5, "%E %E %E %E %E %E %E %E %E %E\n", gr.tab[jr].r, newtm->sig_ster[jr], newtm->sig_cor[jr],
                    newtm->mz[jr], newtm->hd_ster[jr],
                    newtm->Ptot[jr], newtm->Pgas[jr], newtm->nu_ster[jr],
                    newtm->vr_dysk[jr], newtm->vr_cor[jr]);
        } else {

	  fprintf(plik5, "%E %E %E %E %E %E %E %E\n", gr.tab[jr].r, newtm->sig_ster[jr], newtm->hd_ster[jr], newtm->Ptot[jr], newtm->Pgas[jr], newtm->nu_ster[jr], newtm->vr_dysk[jr], newtm->Alpha22[jr]);
        }

    }

#endif

    if(rank == 0)
        fclose(plik5);

    return 1;

}


#ifdef HDF5_USED
/**************************************************************/
int WriteCoordsHDF(const char *FileName)
/**************************************************************/
{
    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[2];                 /* dataset dimensions */
    double     *data;                    /* pointer to data buffer to write */
    hsize_t		count[2];	          /* hyperslab selection parameters */
    hsize_t		offset[2];
    int         i;
    herr_t		status;
	char 		fullFileName[256];

	if(rank == 0)
		printf("Writing HDF file %s\n", FileName);

	sprintf(fullFileName, "hdf/%s", FileName);
	mkdir("hdf", 0777);

    file_id = H5Fcreate(fullFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(file_id < 0) printf("HDF Error 2\n");
   
	//All
	{
		dimsf[0] = PROZ * ZROZ;
		dimsf[1] = 2;
		filespace = H5Screate_simple(2, dimsf, NULL); 
		/*
		* Create the dataset with default properties and close filespace.
		*/
		dset_id = H5Dcreate(file_id, "/coordsXY", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(filespace);

		/*
		 * Initialize data buffer 
		 */
		data = (double *) malloc(dimsf[0] * dimsf[1] * sizeof(double));

		i = 0;
		for (int jr = 0; jr < PROZ; jr++)
		  for (int z = 0; z < ZROZ; z++)
		  {
			data[i++] = log10(gr.tab[jr].r*Rschw);
			if(gr.z[z] > 0.0)
				data[i++] = log10(gr.z[z]);
			else if(gr.z[z] < 0.0)	
				data[i++] = -log10(-gr.z[z]);
			else	
				data[i++] = 0.0;
		  }

		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		free(data);	

		/*
		 * Close/release resources.
		 */
		H5Dclose(dset_id);
	}
	
	
	//r
	{
		dimsf[0] = PROZ;
		filespace = H5Screate_simple(1, dimsf, NULL); 
		/*
		* Create the dataset with default properties and close filespace.
		*/
	//	strcpy(FieldName, "mesh");
		hid_t par_id = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(par_id, 1, dimsf);	
		
		dset_id = H5Dcreate(file_id, "/r", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, par_id, H5P_DEFAULT);
		H5Pclose(par_id);
		H5Sclose(filespace);

		/*
		 * Initialize data buffer 
		 */
		data = (double *) malloc(PROZ * sizeof(double));

		for (int jr = 0; jr < PROZ; jr++)
			data[jr] = gr.tab[jr].r*Rschw;

		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		free(data);	

		/*
		 * Close/release resources.
		 */
		H5Dclose(dset_id);
	}

	//phi
#if 0
	{
		dimsf[0] = konf.phi_size;
		filespace = H5Screate_simple(1, dimsf, NULL); 
		/*
		* Create the dataset with default properties and close filespace.
		*/
		dset_id = H5Dcreate(file_id, "/phi", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Pclose(par_id);
		H5Sclose(filespace);

		/*
		 * Initialize data buffer 
		 */
		data = (double *) malloc(konf.phi_size * sizeof(double));

		for (int p = 0; p < konf.phi_size; p++)
		{
			data[p] = 2.0*Pi*2.0/3.0*((double)p/konf.phi_size + Pi/3.0);
			printf("p=%02d dp=%E\n", p, data[p]);
		}	

		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		free(data);	

		/*
		 * Close/release resources.
		 */
		H5Dclose(dset_id);
	}
#endif	

	//z
	{
		dimsf[0] = ZROZ;
		filespace = H5Screate_simple(1, dimsf, NULL); 
		/*
		* Create the dataset with default properties and close filespace.
		*/
	//	strcpy(FieldName, "mesh");
		hid_t par_id = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(par_id, 1, dimsf);	
		
		dset_id = H5Dcreate(file_id, "/z", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, par_id, H5P_DEFAULT);
		H5Pclose(par_id);
		H5Sclose(filespace);

		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, gr.z);
		/*
		 * Close/release resources.
		 */
		H5Dclose(dset_id);
	}

	H5Fclose(file_id);	
    return 1;
}


/**************************************************************/
int WriteOutputHDF(double time, int ntime, SOLINTIME *newtm)
/**************************************************************/
{
    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[2];                 /* dataset dimensions */
    double     *data;                    /* pointer to data buffer to write */
    hsize_t		count[2];	          /* hyperslab selection parameters */
    hsize_t		offset[2];
    hid_t		plist_id;                 /* property list identifier */
    int         i;
    herr_t		status;
    char 		FileName[256];
    char 		fullFileName[260];

    /* 
     * Set up file access property list with parallel I/O access
     */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
	 if(plist_id < 0) printf("HDF Error 1\n");
	 
#ifdef MPI_USED
     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
    /*
     * Create a new file collectively and release property list identifier.
     */
    sprintf(FileName, "hdf%05d.h5", ++konf.lasthdf);
	sprintf(fullFileName, "hdf/%s", FileName);

	if(rank == 0 && konf.lasthdf == 1)
		WriteCoordsHDF("coords.h5");

	if(rank == 0)
		printf("Writing HDF file %s\n", FileName);
	
    file_id = H5Fcreate(fullFileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	if(file_id < 0) printf("HDF Error 2\n");

    H5Pclose(plist_id);
   
	if(rank == 0)
		write_xdmf_xml(time, konf.lasthdf, FileName);

    /*
     * Create the dataspace for the dataset.
     */
    int minjr = (konf.ji == 1 ? 0 : konf.ji);
    int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

	{
		dimsf[0] = PROZ;
		dimsf[1] = ZROZ;
		filespace = H5Screate_simple(2, dimsf, NULL); 

		/*
		* Create the dataset with default properties and close filespace.
		*/
	//	strcpy(FieldName, "mesh");
//		hid_t par_id = H5Pcreate(H5P_DATASET_CREATE);
//		H5Pset_chunk(par_id, 2, dimsf);	
//		status = H5Pset_deflate (par_id, 6);
		
		dset_id = H5Dcreate(file_id, "/density", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT /*par_id*/, H5P_DEFAULT);
//		H5Pclose(par_id);
		H5Sclose(filespace);

		/* 
		 * Each process defines dataset in memory and writes it to the hyperslab
		 * in the file.
		 */
		count[0] = (maxjr - minjr + 1);
		count[1] =  ZROZ;
		
		offset[0] = minjr;
		offset[1] = 0;
		memspace = H5Screate_simple(2, count, NULL);

		/*
		 * Select hyperslab in the file.
		 */
		filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

	/*
		 * Initialize data buffer 
		 */
		data = (double *) malloc(count[0] * count[1] * sizeof(double));

		i = 0;
		const double SMALLZ = 1.0e-6;
//		 for (int p = 0; p < konf.phi_size; p++)
		for (int jr = minjr; jr <= maxjr; jr++)
		  for (int z = 0; z < ZROZ; z++)
		  {
			if(newtm->hd_ster[jr] >= fabs(gr.z[z]) && newtm->hd_ster[jr] > SMALLZ)
				data[i++] = log10(newtm->sig_ster[jr]/newtm->hd_ster[jr]);
	//	data[i++] = newtm->sig_ster[jr]/newtm->hd_ster[jr];
			else	
				data[i++] = 0.0;
//			data[i++] = 1.0e-20;
		  }

#ifdef MPI_USED
		/*
		 * Create property list for collective dataset write.
		 */
		plist_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
				  plist_id, data);
		H5Pclose(plist_id);
#else
	
		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data);
//		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#endif
		free(data);	

		/*
		 * Close/release resources.
		 */
		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Sclose(memspace);
	}


	{
		dimsf[0] = PROZ;
		dimsf[1] = ZROZ;
		filespace = H5Screate_simple(2, dimsf, NULL); 

		/*
		* Create the dataset with default properties and close filespace.
		*/
	//	strcpy(FieldName, "mesh");
//		hid_t par_id = H5Pcreate(H5P_DATASET_CREATE);
//		H5Pset_chunk(par_id, 2, dimsf);	
//		status = H5Pset_deflate (par_id, 6);
		
		dset_id = H5Dcreate(file_id, "/temperature", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT /*par_id*/, H5P_DEFAULT);
//		H5Pclose(par_id);
		H5Sclose(filespace);

		/* 
		 * Each process defines dataset in memory and writes it to the hyperslab
		 * in the file.
		 */
		count[0] = (maxjr - minjr + 1);
		count[1] =  ZROZ;
		
		offset[0] = minjr;
		offset[1] = 0;
		memspace = H5Screate_simple(2, count, NULL);

		/*
		 * Select hyperslab in the file.
		 */
		filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

	/*
		 * Initialize data buffer 
		 */
		data = (double *) malloc(count[0] * count[1] * sizeof(double));

		i = 0;
		const double SMALLZ = 1.0e-6;
//		 for (int p = 0; p < konf.phi_size; p++)
		for (int jr = minjr; jr <= maxjr; jr++)
		  for (int z = 0; z < ZROZ; z++)
		  {
			if(newtm->hd_ster[jr] >= fabs(gr.z[z]) && newtm->hd_ster[jr] > SMALLZ)
				data[i++] = log10(newtm->temp_ster[jr]);

	//		data[i++] = newtm->temp_ster[jr];
			else	
				data[i++] = 0.0;
//			data[i++] = 1.0e-20;
		  }

#ifdef MPI_USED
		/*
		 * Create property list for collective dataset write.
		 */
		plist_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
				  plist_id, data);
		H5Pclose(plist_id);
#else
	
		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data);
//		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#endif
		free(data);	

		/*
		 * Close/release resources.
		 */
		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Sclose(memspace);
	}

	H5Fclose(file_id);	
    return 1;
}

void write_xdmf_xml(double time, int hdfNo, const char *hdf_file)
{
    FILE *xmf = 0;

    /*
     * Open the file and write the XML description of the mesh..
     */
    xmf = fopen("hdf/DZWM_model.xmf", hdfNo == 1 ? "w" : "r+");

	if(hdfNo == 1)
	{
		fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
		fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
		fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
		fprintf(xmf, " <Domain>\n");
		fprintf(xmf, " <Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
	}
	else
	{
		fseek(xmf, -strlen(" </Grid>\n" " </Domain>\n" "</Xdmf>\n"), SEEK_END);
	}
	
    fprintf(xmf, "   <Grid Name=\"Model for t=%.1lf\" GridType=\"Uniform\">\n", time);
    fprintf(xmf, "   <Time Value=\"%.1lf\"/>\n", time);
//    fprintf(xmf, "   <Information Name=\"Dims\" Value=\"%d\"/>\n", 2);
//    fprintf(xmf, "   <Information Name=\"Dim0\" Value=\"%s\"/>\n", "log10(r)");
//    fprintf(xmf, "   <Information Name=\"Unit0\" Value=\"%s\"/>\n", "m");
//    fprintf(xmf, "   <Information Name=\"Dim1\" Value=\"%s\"/>\n", "log10(z)");
//    fprintf(xmf, "   <Information Name=\"Unit1\" Value=\"%s\"/>\n", "m");	
    fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", PROZ, ZROZ);
    fprintf(xmf, "     <Geometry GeometryType=\"XY\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", PROZ*ZROZ, 2);
    fprintf(xmf, "        coords.h5:/coordsXY\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", PROZ, ZROZ);
    fprintf(xmf, "        %s:/density\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"Temperature\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", PROZ, ZROZ);
    fprintf(xmf, "        %s:/temperature\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");


    fprintf(xmf, " </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}
#endif //#ifdef HDF5_USED
 

/**************************************************************/
int FinalizeOutput(void)
/**************************************************************/
{
    int minjr = (konf.ji == 1 ? 0 : konf.ji);
    int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

    if(konf.lsig)
        for (int jr = minjr; jr <= maxjr; jr++) {
            fclose(sigtefiles[jr]);
            //        fclose(mdoty[jr]);
        }

    if(konf.llum && rank == 0)
        fclose(lumfile);

    return 1;
}

/**************************************************************/
int WriteKrzyweS(void)
/**************************************************************/
{
    register int jr;
    register int ii;
    int minjr = (konf.ji == 1 ? 0 : konf.ji);
    int maxjr = (konf.jo == PROZ - 2 ? PROZ-1 : konf.jo);

    char s_file[256];

    FILE *plik;

    for (jr = minjr; jr <= maxjr; jr++) {
        //drukowanie krzywych stacjonarnych
        sprintf(s_file, "sigma.%03d", jr);

        if((plik=fopen(s_file, "w"))==NULL) {
            printf("Cannot open file %s\n", s_file);
            return 0;
        }

        //	       printf("Printing S-curve for r=%lf; file = %s\n", gr.tab[jr].r, s_file);

        for(ii= 0; ii < MROZ; ii++) {
            fprintf(plik,"%E %E %E %E\n", gr.tab[jr].sigma[ii], gr.tab[jr].Tef[ii],
                    gr.tab[jr].temper[ii], gr.tab[jr].Mdotstab[ii]);
        }
        fclose(plik);

        //    printf("jr=%d\n", jr);

    }

    //    printf("jr=%d\n", jr);



    if(rank == 0)
        printf("S-curves printed to files; last file for r=%lf\n", gr.tab[PROZ-1].r);

    return 1;
}

/**************************************************************/
int WriteOneKrzywaS(int jr)
/**************************************************************/
{
    register int ii;

    char s_file[256];

    FILE *plik;

    //drukowanie krzywych stacjonarnych
    sprintf(s_file, "sigma.%03d", jr);

    if((plik=fopen(s_file, "w"))==NULL) {
        printf("Cannot open file %s\n", s_file);
        return 0;
    }

    //	       printf("Printing S-curve for r=%lf; file = %s\n", gr.tab[jr].r, s_file);

    for(ii= 0; ii < MROZ; ii++) {
        fprintf(plik,"%E %E %E %E\n", gr.tab[jr].sigma[ii], gr.tab[jr].Tef[ii],
                gr.tab[jr].temper[ii], gr.tab[jr].Mdotstab[ii]);
    }
    fclose(plik);

    //    printf("jr=%d\n", jr);

    printf("S-curve printed to file %s for prom = %E\n", s_file, gr.tab[jr].r);

    return 1;
}


/**************************************************************/
int RestartDump(SOLINTIME *curtm, double *dt)
/**************************************************************/
{
    FILE *bf;

    char dump_name[256];

    sprintf(dump_name, "dump_%d_%d.dat.%06d", rank, size, ++konf.nrestore);

    if(rank == 0)
        printf("Restart dump no %d at time %E\n", konf.nrestore, ptime);

    if((bf = fopen(dump_name,"wb"))==NULL) {
        printf("Cannot open file %s\n", dump_name);
        return 0;
    }

    SaveFilePositions();           //Zapamietujemy aktualne pozycje w plikach

    //Zapisujemy akualne dane

    if(fwrite(&ptime, sizeof(ptime), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }

    if(fwrite(&ntime, sizeof(ntime), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }

    if(fwrite(dt, sizeof(double), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }

    if(fwrite(&konf, sizeof(KONFIG), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }


    if(fwrite(&gr, sizeof(GRID), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }

    if(fwrite(curtm, sizeof(SOLINTIME), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }

    if(fwrite(&solve, sizeof(SOLVESTRUCT), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }

    if(fwrite(&filepos, sizeof(filepos), 1, bf) != 1) {
        printf("Error writing file %s\n", dump_name);
        return 0;
    }


    fclose(bf);

    if(rank == 0)
        printf("Restart dump saved\n");

    return 1;
}


/**************************************************************/
int RestartRestore(int dump_num, SOLINTIME *curtm, double *dt)
/**************************************************************/
{
    FILE *bf;

    char dump_name[256];

    sprintf(dump_name, "dump_%d_%d.dat.%06d", rank, size, dump_num);

    if((bf = fopen(dump_name,"rb"))==NULL) {
        printf("Cannot open file %s\n", dump_name);
        return 0;
    }


    //Zapisujemy dane z dumpa

    if(fread(&ptime, sizeof(ptime), 1, bf) != 1) {
        printf("Error reading file%s\n", dump_name);
        return 0;
    }

    if(fread(&ntime, sizeof(ntime), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }

    if(fread(dt, sizeof(double), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }

    if(fread(&konf, sizeof(KONFIG), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }


    if(fread(&gr, sizeof(GRID), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }

    if(fread(curtm, sizeof(SOLINTIME), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }

    if(fread(&solve, sizeof(SOLVESTRUCT), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }

    if(fread(&filepos, sizeof(filepos), 1, bf) != 1) {
        printf("Error reading file %s\n", dump_name);
        return 0;
    }

    fclose(bf);

    Mext   = konf.Mdot*Msun/Year;
    Rschw  = 2.0*Graw*konf.Mass/Crad/Crad;
    Ledd   = 1.25e39*(konf.Mass/10.0/Msun);

    return 1;
}

/**************************************************************/
void dataio(double ptime, int ntime, double dt, SOLINTIME *curtm)
/**************************************************************/
{


    if (konf.llum > 0 && (ptime  >= (konf.tlum + konf.dtlum) || (konf.deltalum > 0.0 && ptime  >= konf.tlum))) {

        double lum, lum3, lumcor;

        lum = curtm->lum2;
        lum3 = curtm->lum3;
        lumcor = curtm->lumcor;

#ifdef MPI_USED
        //Pobieramy sume jasnosci
        if (size > 1) {
            double isum[3], osum[3];
            isum[0] = lum;
            isum[1] = lum3;
            isum[2] = lumcor;

            MPI_Reduce( isum, osum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

            lum    = osum[0];
            lum3   = osum[1];
            lumcor = osum[2];
        }
#endif
        lum = log10(lum);
        lum3 = log10(lum3);

        if(konf.wymiana_masy)
            lumcor = log10(lumcor);

        if (ptime  >= (konf.tlum + konf.dtlum)) {
            filepos.lastlum = lum;
            filepos.lastlumcor = lumcor;

            if(rank == 0)
                WriteOutputLum(ptime, ntime, lum, lum3, lumcor);
            konf.tlum += konf.dtlum;
        } else if(konf.deltalum > 0.0 && (fabs(lum - filepos.lastlum) >= konf.deltalum || fabs(lumcor - filepos.lastlumcor) >= konf.deltalum)) {

            filepos.lastlum = lum;
            filepos.lastlumcor = lumcor;

            if(rank == 0)
                WriteOutputLum(ptime, ntime, lum, lum3, lumcor);
        }
    }

    if (konf.lsig > 0 && ptime  >= (konf.tsig + konf.dtsig)) {

        WriteOutputSigte(ptime, ntime, curtm);
        konf.tsig += konf.dtsig;
    }

    if (konf.lrad > 0 && ptime  >= (konf.trad + konf.dtrad) && ptime <= (konf.tradend + dt)) {
        WriteOutputRad(ptime, ntime, curtm);
        konf.trad += konf.dtrad;
    }

#ifdef HDF5_USED
    if (konf.lhdf > 0 && ptime  >= (konf.thdf + konf.dthdf)) {
        WriteOutputHDF(ptime, ntime, curtm);
        konf.thdf += konf.dthdf;
    }
#endif    

    if (konf.lrestore > 0 && time(NULL) >= konf.trestore + konf.dtrestore) {
        konf.trestore += konf.dtrestore;
        RestartDump(curtm, &dt);
    }


}

