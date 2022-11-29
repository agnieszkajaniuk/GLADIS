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
#include "mdothold.h"
#include "spline.h"


/**************************************************************/
MdotHold :: MdotHold(const char *filename)
/**************************************************************/
{
    tm  = NULL;
    md  = NULL;
    count = 0;

    readFile(filename);
}

/**************************************************************/
MdotHold :: ~MdotHold()
/**************************************************************/
{
    delete[] tm;
    delete[] md;

}

/**************************************************************/
double MdotHold :: findMdot(double time_)
/**************************************************************/
{
    double ret = 0.0;
    if(!count)
        return 0.0;

    if(time_ <= tm[0])
        ret = md[0];
    else if(time_ >= tm[count-1])
        ret = md[count-1];
    else
        ret = splineInterp(tm, md, count, time_);

    printf("findMdot: ret val = %E for time %E\n", ret, time_);
    return ret;
}


/**************************************************************/
void MdotHold :: readFile(const char *filename)
/**************************************************************/
{
    FILE *fp;

    int   cnt = 0;
    char  buf[256];
    int   tmp;


    if((fp = fopen(filename, "rt")) == NULL)
    {

        printf("Can not open file %s\n", filename);
        return;
    }

    //count number of lines

    while(!feof(fp))
        if(fgets(buf, sizeof(buf), fp))
            cnt++;


    rewind(fp);

    tm = new double[cnt];
    md = new double[cnt];

    for(int i = 0; i < cnt; i++)
    {
        fgets(buf, sizeof(buf), fp);

        sscanf(buf, "%d %lf %lf", &tmp, &tm[i], &md[i]);

        //printf("Read: %s\n", buf);
        //printf("Scan: %E %lf\n", tm[i], md[i]);

        tm[i] = tm[i] * 3.157E+7;
        //    printf("z pliku: md= %E\n", md[i]);

        md[i] = pow(10.0, md[i])/Ledd;

        //  printf("po przeliczeniu: md= %E\n", md[i]);
        printf("MdotHold[%d] = time: %E val: %E\n", i, tm[i], md[i]);

        //   getc(stdin);
    }

    count = cnt;


    fclose(fp);
}

/*****************************************************************/
double splineInterp(double *tabx, double *taby, int count, double valat)
/*****************************************************************/
{
    double ret = 0.0;

    spline_overhauser_val(1, count, tabx, taby, valat, &ret);

    return ret;
}
