
#if !defined(MDOTHOLD_H__INCLUDED_)
#define MDOTHOLD_H__INCLUDED_


//Tabularized mdot
class MdotHold {

private:
    double *tm;       //time
    double *md;       //mdot

    int count;        //number of points

    void readFile(const char *filename);

public:
    MdotHold(const char *filename);
    ~MdotHold();
    double findMdot(double time_);
};


double splineInterp(double *tabx, double *taby, int count, double valat);

#endif
