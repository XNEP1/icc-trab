#ifndef __METODOS_H__
#define __METODOS_H__

enum condicionador{
    PRECONDIONADOR_JACOBI,
    PRECONDIONADOR_NENHUM
};

#define ZERO_EPSILON 1e-15

int conjGradient(SisDiag_t *SD, real_t *x, real_t erro, int maxit, int tipoCondicionador, rtime_t *tempoInter, rtime_t *tempoPC, real_t *maxNormaInter);

real_t normaL2Residuo(SisDiag_t *SD, real_t *x);
void calResiduo(SisDiag_t *SD, real_t *x, real_t *r);

real_t normaL2ResiduoSDSimetrico(SisDiag_t *SDSimetrico, real_t *x);
void calResiduoSDSimetrico(SisDiag_t *SDSimetrico, real_t *x, real_t *r);

#endif // __METODOS_H__