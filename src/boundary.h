
#if !defined(BOUNDARY_H__INCLUDED_)
#define BOUNDARY_H__INCLUDED_

void InitOuterBoundary(void);
void DeInitOuterBoundary(void);

void Boundary(SOLINTIME *curtm);

int OuterBoundary(double *mpar, SOLINTIME *curtm, double time_);

double CalcUrandAlfa(int jr, double urandalfa, double time_, RANDOMDEF *rd, SOLINTIME *newtm);
double Flicker(int jr, double urandalfa, double prom, double time_, double H);

#endif
