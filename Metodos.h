#ifndef __METODOS_H__
#define __METODOS_H__

int conjGradient(SisDiag_t *SD, real_t *x, real_t *M, real_t erro, int maxit, rtime_t *tempoInter, real_t *maxNormaInter);

real_t normaL2Residuo(SisDiag_t *SD, real_t *x);
void calResiduo(SisDiag_t *SD, real_t *x, real_t *r);

real_t normaL2ResiduoSDSimetrico(SisDiag_t *SDSimetrico, real_t *x);
void calResiduoSDSimetrico(SisDiag_t *SDSimetrico, real_t *x, real_t *r);

#endif // __METODOS_H__