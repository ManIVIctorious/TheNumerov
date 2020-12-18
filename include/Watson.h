#ifndef _WATSON_H
#define _WATSON_H

void   init_watson(settings* prefs);
double exec_watson(double*** mu, double** zeta, int* nq, double dq, double** q, int index, int* shift);
void   free_watson(void);

#endif
