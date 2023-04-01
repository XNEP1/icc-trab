/*
    Trabalho 1 de ICC
    Por Bruno Krügel
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "utils.h"
#include "sisDiag.h"
#include "Metodos.h"

#define UF 4 // fator de unrolling
#define BK 8 // Tamanho dos blocos

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void transposta(SisDiag_t *sl, SisDiag_t *sl_Transposta){
    int n = sl->n;
    int q = (sl->k - 1)/2;
    int j_start;

    for(int i=0; i< n; i++){

        j_start = i-q;
        j_start = IS_POSITIVE(j_start) * j_start;

        for(int j=j_start; j < n && j<=i+q; j++)
            sl_Transposta->A[i][j] = sl->A[j][i];
    }
}

void mmm_diag_blocking(SisDiag_t * restrict A, SisDiag_t * restrict B, SisDiag_t * restrict R){

    int iend, jend, kend, ii, jj, kk, jstart;
    int n = A->n;
    int n_bk = n/BK;
    int l = n - (n % BK);
    double r;
    int q_start;
    int q = (R->k - 1)/2;

    int ki_start, kj_start;
    int d = (A->k - 1)/2;

    for(int ii=0; ii<n; ii+=BK){
        iend = MIN(ii+BK,n);

        for(jj=0; jj<n; jj+=BK){
            jend = MIN(jj+BK,n);

            for(kk=0; kk<n; kk+=BK){
                kend = MIN(kk+BK, n);

                for(int i=ii; i<iend; i++){

                    q_start = i - q;
                    q_start = IS_POSITIVE(q_start) * q_start;
                    jstart = MAX(jj, q_start);

                    kj_start = i-d;
                    kj_start = IS_POSITIVE(kj_start) * kj_start;
                    kj_start = MAX(kk, kj_start);

                    for(int j=jstart; j<jend && j<=i+q; j++){
                        r = 0.0;

                        // Serviria pra tirar o if(IS_ZERO(B, k, j))
                        // Mas acaba que fica mais demorado com isso.
                        ki_start = j-d;
                        ki_start = IS_POSITIVE(ki_start) * ki_start;
                        ki_start = MAX(kj_start, ki_start);

                        for(int k=ki_start; k<kend && k<=i+d && k<=j+d; k++){
                        // for(int k=kj_start; k<kend && k<=i+d && k<=j+d; k++){

                            // if(IS_ZERO(B, k, j))
                                // continue;

                            r +=  A->A[i][k] * B->A[k][j];
                        }

                        R->A[i][j] += r;
                    }
                }
            }
        }
    }
}

void mmm_diag(SisDiag_t * restrict A, SisDiag_t * restrict B, SisDiag_t * restrict R){

    int iend, jend, kend, ii, jj, kk, jstart;
    int n = A->n;
    int n_bk = n/BK;
    int l = n - (n % BK);
    double r;
    int j_start;
    int q = (R->k - 1)/2;

    int ki_start, kj_start;
    int d = (A->k - 1)/2;

    for(int i=0; i<n; i++){

        j_start = i - q;
        j_start = IS_POSITIVE(j_start) * j_start;

        kj_start = i-d;
        kj_start = IS_POSITIVE(kj_start) * kj_start;

        for(int j=j_start; j<n && j<=i+q; j++){
            r = 0.0;

            ki_start = j-d;
            ki_start = IS_POSITIVE(ki_start) * ki_start;
            ki_start = MAX(kj_start, ki_start);

            for(int k=ki_start; k<n && k<=i+d && k<=j+d; k++){

                r +=  A->A[i][k] * B->A[k][j];
            }

            R->A[i][j] += r;
        }
    }
}

void mmv_diag_unroll(SisDiag_t * restrict A, real_t * restrict B, real_t * restrict R){
    int iend, jend, kend, ii, jj, kk, jstart;
    int q_start;
    int n = A->n;
    int d = (A->k - 1)/2;
    double temp1, temp2, temp3, temp4;

    int m = n - (n % UF);
    for(int i=0; i<m; i+=UF){
        temp1 = 0.0;
        temp2 = 0.0;
        temp3 = 0.0;
        temp4 = 0.0;
        for(int j=0; j<n; j++){
            temp1 += GET_A(A, i, j)   * B[j];
            temp2 += GET_A(A, i+1, j) * B[j];
            temp3 += GET_A(A, i+2, j) * B[j];
            temp4 += GET_A(A, i+3, j) * B[j];
        }
        R[i]   = temp1;
        R[i+1] = temp2;
        R[i+2] = temp3;
        R[i+3] = temp4;
    }
    // Resto
    for(int i=m; i<n; i++){
        q_start = i-d;
        q_start = IS_POSITIVE(q_start) * q_start;
        temp1 = 0.0;
        for(int j=q_start; j < n && j<=i+d; j++)
            temp1 += A->A[i][j] * B[j];
        R[i] = temp1;
    }
}

void mmv_diag(SisDiag_t * restrict A, real_t * restrict B, real_t * restrict R){
    int d = (A->k - 1)/2;
    int j_start;
    double temp;

    for(int i=0; i < A->n; i++){

        j_start = i-d;
        j_start = IS_POSITIVE(j_start) * j_start;

        temp = 0.0;
        for(int j=j_start; j < A->n && j<=i+d; j++)
            temp += A->A[i][j] * B[j];
        R[i] = temp;
    }
}

/* 
    A^T * Ax = A^T * b 
*/
void simetrizarSL(SisDiag_t *sl, SisDiag_t *slSimetrico){

    double r;
    int q = (slSimetrico->k - 1)/2;

    int n = sl->n;


    SisDiag_t *sl_Transposta = alocaSisDiag(sl->n, sl->k);

    transposta(sl, sl_Transposta);

    /* Multiplicando A^T * A */
    #ifdef UNROLL
    mmm_diag_blocking(sl_Transposta, sl, slSimetrico);
    #else
    mmm_diag(sl_Transposta, sl, slSimetrico);
    #endif

    /* Multiplicando A^T * b  */
    #ifdef UNROLL
    mmv_diag_unroll(sl_Transposta, sl->b , slSimetrico->b );
    #else
    mmv_diag(sl_Transposta, sl->b , slSimetrico->b );
    #endif

    liberaSisDiag(sl_Transposta);
}

int conjGradient(SisDiag_t *SD, real_t *x, real_t erro, int maxit, int tipoCondicionador ,rtime_t *tempoInter, rtime_t *tempoPC, real_t *maxNormaInter){
    rtime_t tempoTodasInteracoes;
    rtime_t tempoAux;

    SisDiag_t *sdSim = alocaSisDiag(SD->n, 2*SD->k - 1);
    unsigned int n = sdSim->n;
    int q = (sdSim->k - 1) / 2;
    real_t *r = malloc(n * sizeof(real_t));
    real_t *p = malloc(n * sizeof(real_t));
    real_t *A_p = malloc(n * sizeof(real_t));
    real_t *z = malloc(n * sizeof(real_t));
    real_t *precondicionador = (real_t *) malloc(n * sizeof(real_t));
    real_t *x_old = malloc(n * sizeof(real_t));
    unsigned int interacoes = 0;
    real_t e = erro * 2;
    real_t a = 0.0;
    real_t b = 0.0;
    real_t temp1, temp2, temp3;
    int m;

    int i, j, k;

    if(!sdSim || !r || !p || !z || !x_old ){
        return -2;
    }

    tempoAux = timestamp();
    simetrizarSL(SD, sdSim);

    switch (tipoCondicionador){
    case PRECONDIONADOR_JACOBI:
        for(i=0; i<n; i++)
            precondicionador[i] = (1/sdSim->A[i][i]);  // m é a diagonal principal do SL (pré-condicionador de Jacobi)
        break;

    case PRECONDIONADOR_NENHUM:
        for(i=0; i<n; i++)
            precondicionador[i] = 1; // sem pré-condicionador
        break;
    
    default:
        return -3;
        break;
    }
    *tempoPC += timestamp() - tempoAux;

    calResiduo(sdSim, x, r);
    for(i=0; i<n; i++)
        e += r[i] * r[i];
    e = sqrt(e); /* 'e' sempre será positivo. */

    for(i=0; i<n; i++)
        z[i] = precondicionador[i]*r[i];

    for(i=0; i<n; i++)
        p[i] = z[i];

    tempoTodasInteracoes = timestamp();
    while(e > erro && interacoes < maxit ){

        /* Calcula o alfa */
        temp1 = 0.0;
        for(i=0; i<n; i++)
            temp1 += r[i]   * z[i];

        #ifdef UNROLL
        mmv_diag_unroll(sdSim, p, A_p);
        #else
        mmv_diag(sdSim, p, A_p);
        #endif

        temp2 = 0.0;
        for(i=0; i<n; i++){
            temp2 += A_p[i] * p[i];
        }

        if(ABS(temp2) > ZERO_EPSILON){ // sem divisão por 0
            a = temp1/temp2;
        }

        /* Atualiza x */
        for(i=0; i<n; i++){
            x_old[i] = x[i];
            x[i] = x[i] + a*p[i];
        }

        /* Atualiza r */
        for(i=0; i<n; i++){
            r[i] = r[i] - a*A_p[i];
        }

        /* Atualiza z */
        for(i=0; i<n; i++)
            z[i] = precondicionador[i]*r[i];

        /* Calcula beta */
        b = 0.0;
        temp2 = 0.0;
        for(i=0; i<n; i++)
            temp2 += r[i] * z[i];


        if(ABS(temp1) > ZERO_EPSILON){ // sem divisão por 0
            b = temp2/temp1;
        }

        /* Atualiza p */
        for(i=0; i<n; i++)
            p[i] = z[i] + b*p[i];

        
        /* Calcula a norma maxima da interação */
        maxNormaInter[interacoes] = x[0] - x_old[0];
        int maxIndex = 0;
        for(i=1; i<n; i++){
            temp1 = fabs(x[i] - x_old[i]); 
            if(maxNormaInter[interacoes] < temp1){
                maxNormaInter[interacoes] = temp1;
                maxIndex = i;
            }
        }


        /* criterio de parada */
        e = maxNormaInter[interacoes] / fabs(x[maxIndex]);

        interacoes = interacoes + 1;
    }

    tempoTodasInteracoes = timestamp() - tempoTodasInteracoes;
    free(r);
    free(z);
    free(p);
    free(A_p);
    free(x_old);
    liberaSisDiag(sdSim);
    *tempoInter = tempoTodasInteracoes / interacoes;
    return interacoes;
}

real_t normaL2Residuo(SisDiag_t *SD, real_t *x){
    int i;
    real_t *r = malloc(SD->n * sizeof(real_t));
    real_t e = 0.0;

    calResiduo(SD, x, r);

    for(i=0; i<SD->n; i++)
        e += r[i] * r[i];
    e = sqrt(e); /* 'e' sempre será positivo. */

    free(r);
    return e;
}

void calResiduo(SisDiag_t *SD, real_t *x, real_t *r){
    int i, j;
    int j_start;

    int q = (SD->k - 1) / 2;
    for(int i=0; i < SD->n; i++){

        j_start = i - q;
        j_start = IS_POSITIVE(j_start) * j_start;

        r[i] = 0.0;
        for(int j=j_start; j < SD->n && j<=i+q; j++) // Só percorre j's que não são de diagonais nulas
            r[i] -= SD->A[i][j] * x[j];

        r[i] += SD->b[i];
    }
}