#ifndef _WATSON_H
#define _WATSON_H

void   init_watson(settings* prefs, double dq);
double exec_watson_2d(double*** mu, double** zeta, double** q, int index, int* shift);
double exec_watson_3d(double*** mu, double** zeta, double** q, int index, int* shift);
void   free_watson(void);

#endif
