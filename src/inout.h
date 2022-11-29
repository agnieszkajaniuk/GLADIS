
#if !defined(INOUT_H__INCLUDED_)
#define INOUT_H__INCLUDED_


int InitializeOutput(int restart);
int WriteOutputLum(double time, int ntime, double lum2, double lum3, double lumcor);
int WriteOutputSigte(double time, int ntime, SOLINTIME *newtm);
int WriteOutputRad(double time, int ntime, SOLINTIME *newtm);
int WriteOutputHDF(double time, int ntime, SOLINTIME *newtm);
int FinalizeOutput(void);

int WriteKrzyweS(void);
int WriteOneKrzywaS(int jr);


int RestartDump(SOLINTIME *curtm, double *dt);
int RestartRestore(int dump_num, SOLINTIME *curtm, double *dt);

void dataio(double ptime, int ntime, double dt, SOLINTIME *curtm);

#endif
