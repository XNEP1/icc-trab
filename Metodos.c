/*
    Trabalho 1 de ICC
    Por Bruno Krügel
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sisDiag.h"
#include "Metodos.h"

int conjGradient(SisDiag_t *SD, real_t *x, real_t *M, real_t erro, int maxit, rtime_t *tempoInter, real_t *maxNormaInter){
    rtime_t tempoTodasInteracoes;

    unsigned int n = SD->n;
    real_t *r = malloc(n * sizeof(real_t));
    real_t *p = malloc(n * sizeof(real_t));
    real_t *z = malloc(n * sizeof(real_t));;
    real_t *x_old = malloc(n * sizeof(real_t));
    unsigned int interacoes = 0;
    real_t e = erro * 2;
    real_t a;
    real_t b;
    real_t q1, q2, q3;

    int i, j, i_1, j_1, k;

    calResiduoSDSimetrico(SD, x, r);
    for(i=0; i<n; i++)
        e += r[i] * r[i];
    e = sqrt(e); /* 'e' sempre será positivo. */

    for(i=0; i<n; i++)
        z[i] = (1/M[i])*r[i];

    for(i=0; i<n; i++)
        p[i] = z[i];

    tempoTodasInteracoes = timestamp();
    while(e > erro && interacoes < maxit ){

        /* Calcula o alfa */
        q1 = 0.0;
        for(i=0; i<n; i++)
            q1 += r[i] * z[i];

        q2 = 0.0;
        for(j=0; j<n; j++){
            q3 = 0.0;
            for(i=0; i<n; i++){

                k = j - i + (SD->k - 1);  // Coluna onde o numero está
                if(k >= (SD->k*2 - 1) || k < 0) // é zero.
                    continue;

                if(k <= (SD->k - 2)){     // Está em uma diagonal inferior, então os indices precisam ser ajustado para sua equivalencia nas diagonais superiores.
                    i_1 = i + k - (SD->k - 1);
                    j_1 = SD->k - k - 1;

                    q3 += (p[i] * SD->A[i_1][j_1]);
                }
                else{   // O número está em uma das diagonais superiores.
                    j_1 = k - (SD->k - 1);
                    q3 += (p[i] * SD->A[i][j_1]);
                }
            }
            q2 += q3 * p[j];
        }

        a = q1/q2;


        /* Atualiza x */
        for(i=0; i<n; i++){
            x_old[i] = x[i];
            x[i] = x[i] + a*p[i];
        }

        /* Atualiza r */
        for(j=0; j<n; j++){
            q2 = 0.0;
            for(i=0; i<n; i++){

                k = j - i + (SD->k - 1);  // Coluna onde o numero está
                if(k >= (SD->k*2 - 1) || k < 0) // é zero.
                    continue;

                if(k <= (SD->k - 2)){     // Está em uma diagonal inferior, então os indices precisam ser ajustado para sua equivalencia nas diagonais superiores.
                    i_1 = i + k - (SD->k - 1);
                    j_1 = SD->k - k - 1;

                    q2 += (p[i] * SD->A[i_1][j_1]);
                }
                else{   // O número está em uma das diagonais superiores.
                    j_1 = k - (SD->k - 1);
                    q2 += (p[i] * SD->A[i][j_1]);
                }

            }

            r[j] = r[j] - a*q2;
        }

        for(i=0; i<n; i++)
            z[i] = (1/M[i])*r[i];

        /* Calcula beta */
        b = 0.0;
        q2 = 0.0;
        for(i=0; i<n; i++)
            q2 += r[i] * z[i];
        
        b = q2/q1;

        /* Atualiza p */
        for(i=0; i<n; i++)
            p[i] = z[i] + b*p[i];

        
        /* Calcula a norma maxima da interação */
        maxNormaInter[interacoes] = x[0] - x_old[0];
        for(i=1; i<n; i++){
            q1 = x[i] - x_old[i]; 
            if(maxNormaInter[interacoes] < q1)
                maxNormaInter[interacoes] = q1;
        }


        /* criterio de parada */
        e = maxNormaInter[interacoes];

        interacoes = interacoes + 1;
    }

    tempoTodasInteracoes = timestamp() - tempoTodasInteracoes;

    free(r);
    free(z);
    free(p);
    free(x_old);
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
    int i, j, k;

    for(i=0; i<SD->n; i++){
        r[i] = 0.0;
        for(j=0; j<SD->n; j++){
            k = j - i + (SD->k/2);
            if(k >= SD->k || k < 0) // é zero.
                continue;
            r[i] -= SD->A[i][k] * x[j];
        }
        r[i] += SD->b[i];
    }
}

real_t normaL2ResiduoSDSimetrico(SisDiag_t *SDSimetrico, real_t *x){
    int i;
    real_t *r = malloc(SDSimetrico->n * sizeof(real_t));
    real_t e = 0.0;

    calResiduoSDSimetrico(SDSimetrico, x, r);

    for(i=0; i<SDSimetrico->n; i++)
        e += r[i] * r[i];
    e = sqrt(e); /* 'e' sempre será positivo. */

    free(r);
    return e;
}

void calResiduoSDSimetrico(SisDiag_t *SDSimetrico, real_t *x, real_t *r){
    int i, j, i_1, j_1 , k;

    for(i=0; i<SDSimetrico->n; i++){
        r[i] = 0.0;
        for(j=0; j<SDSimetrico->n; j++){

            k = j - i + (SDSimetrico->k - 1);  // Coluna onde o numero está
            if(k >= (SDSimetrico->k*2 - 1) || k < 0) // é zero.
                continue;

            if(k <= (SDSimetrico->k - 2)){     // Está em uma diagonal inferior, então os indices precisam ser ajustado para sua equivalencia nas diagonais superiores.
                i_1 = i + k - (SDSimetrico->k - 1);
                j_1 = SDSimetrico->k - k - 1;

                r[i] -= SDSimetrico->A[i_1][j_1] * x[j];
            }
            else{   // O número está em uma das diagonais superiores.
                j_1 = k - (SDSimetrico->k - 1);
                r[i] -= SDSimetrico->A[i][j_1] * x[j];
            }

        }
        r[i] += SDSimetrico->b[i];
    }
}