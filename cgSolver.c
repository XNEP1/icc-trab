/*
    Trabalho 1 de ICC
    Por Bruno Krügel
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <likwid.h>

#include "utils.h"
#include "sisDiag.h"
#include "Metodos.h"

/* 
    A^T * Ax = A^T * b 
*/
void simetrizarSL(SisDiag_t *sl, SisDiag_t *slSimetrico){
    int i, j, k_1, k_2, q, p;
    real_t sum;
    q = sl->k / 2;

    /* Multiplicando A^T * A */
    for(i=0; i<sl->n; i++){
        for(j=i; j<sl->n; j++){

            k_1 = j - i;
            if( i >= slSimetrico->n || k_1 >= slSimetrico->k )
                continue; /* Fora da banda de diagonais */

            sum = 0.0;
            for(p=0; p<sl->n; p++){
                k_1 = i - p + q;
                k_2 = j - p + q;
                if(k_1 >= sl->k || k_1 < 0) // é zero.
                    continue;
                if(k_2 >= sl->k || k_2 < 0) // é zero.
                    continue;
                sum += sl->A[p][k_1] * sl->A[p][k_2];
            }

            k_1 = j - i;
            slSimetrico->A[i][k_1] = sum;
        }
    }

    /* Multiplicando A^T * b  */
    for(i=0; i<sl->n; i++){
        sum = 0.0;
        for(p=0; p<sl->n; p++){
            k_1 = i - p + q;
            if(k_1 >= sl->k || k_1 < 0) // é zero.
                continue;
            sum += sl->A[p][k_1] * sl->b[p];
        }
        slSimetrico->b[i] = sum;
    }

}

int main(int argc, char **argv){

    srand(20222);
    LIKWID_MARKER_INIT;

    int parameter_n;
    int parameter_k;
    int parameter_p;
    int parameter_i;
    real_t parameter_e = -1.0;
    FILE *parameter_outputFile; 

    SisDiag_t *slOtm, *slOtmSimetrico;
    real_t *x; // resultado
    real_t *m; // precondicionador
    real_t *maxNormaInter; // inter k
    real_t normaResiduo;
    rtime_t tempoPC = 0.0;
    rtime_t tempoInter = 0.0;
    rtime_t tempoResiduo = 0.0;
    rtime_t tempoAux;
    int i, j;
    int interacoes;
    int status = 0;

    for(i=1; i<argc; i+=2){
        if(!strcmp(argv[i], "-n"))
            parameter_n = atoi(argv[i+1]);
        else if(!strcmp(argv[i], "-k"))
            parameter_k = atoi(argv[i+1]);
        else if(!strcmp(argv[i], "-p"))
            parameter_p = atoi(argv[i+1]);
        else if(!strcmp(argv[i], "-i"))
            parameter_i = atoi(argv[i+1]);
        else if(!strcmp(argv[i], "-e"))
            parameter_e = strtod(argv[i+1], NULL);
        else if(!strcmp(argv[i], "-o"))
            parameter_outputFile = fopen(argv[i+1], "w");
        else{
            perror("Algum dos argumentos estão errados! Execute com:\ncgSolver -n <n> -k <k> -p <p> -i <i> -e <ε> -o <arquivo_saida>\n");
            exit(-1);
        }
    }

    if(!(parameter_n>10) || !(parameter_k>1 && parameter_k%2) || !(parameter_p >= 0) || !(parameter_i >= 0) || !(parameter_outputFile != NULL)){
        perror("Algum dos argumentos não segue os requesitos.\n");
        exit(-1);
    }

    slOtm = alocaSisDiag(parameter_n, parameter_k);
    slOtmSimetrico = alocaSisDiag(parameter_n, ((parameter_k > parameter_n) ? parameter_n : parameter_k));
    x = (real_t *) calloc(parameter_n , sizeof(real_t)); // Todos os valores de x começam sendo 0.0
    m = (real_t *) malloc(parameter_n * sizeof(real_t));
    maxNormaInter = (real_t *) calloc(parameter_i, sizeof(real_t));

    if(!(slOtm) || !(slOtmSimetrico) || !(x) || !(maxNormaInter) || !(m)){
        perror("Erro de alocação de memoria.\n");
        exit(-2);
    }

    gerarCoeficientesSD(slOtm);

    tempoAux = timestamp();
    simetrizarSL(slOtm, slOtmSimetrico);
    tempoPC += timestamp() - tempoAux;

    //prnSisDiag(slOtm);
    //prnSisDiag(slOtmSimetrico);

    tempoAux = timestamp();
    if(parameter_p == 0)
        for(i=0; i<parameter_n; i++)
            m[i] = 1; // sem pré-condicionador
    else
        for(i=0; i<parameter_n; i++)
            m[i] = slOtmSimetrico->A[i][0];  // m é a diagonal principal do SL (pré-condicionador de Jacobi)
    tempoPC += timestamp() - tempoAux;

    char * MARKER_NAME_op1 = markerName("conjGradient", parameter_n);
    LIKWID_MARKER_START(MARKER_NAME_op1);
    interacoes = conjGradient(slOtmSimetrico, x, m, parameter_e, parameter_i, &tempoInter, maxNormaInter);
    LIKWID_MARKER_STOP(MARKER_NAME_op1);

    char * MARKER_NAME_op2 = markerName("calResiduo", parameter_n);
    real_t *r = malloc(parameter_n * sizeof(real_t));
    LIKWID_MARKER_START(MARKER_NAME_op2);
    calResiduo(slOtm, x, r);
    LIKWID_MARKER_STOP(MARKER_NAME_op2);

    tempoAux  = timestamp();
    normaResiduo = normaL2ResiduoSDSimetrico(slOtmSimetrico, x);
    tempoResiduo = timestamp() - tempoAux;

    if(normaResiduo >= 1.0){
        fprintf(stderr, "Residuo muito alto: o metodo pode não ter encontrado uma solução.\n");
        status = 1;
    }

    fprintf(parameter_outputFile, "# bk20 Bruno Krügel\n");
    fprintf(parameter_outputFile, "#\n");
    for(i=0; i<interacoes; i++)
        fprintf(parameter_outputFile, "# inter %d: %.15g\n", i+1, maxNormaInter[i]);
    fprintf(parameter_outputFile, "# residuo: %.15g\n", normaResiduo);
    fprintf(parameter_outputFile, "# Tempo PC: %12g\n", tempoPC);
    fprintf(parameter_outputFile, "# Tempo iter: %12g\n", tempoInter);
    fprintf(parameter_outputFile, "# Tempo iter total: %12g\n", tempoInter * interacoes);
    fprintf(parameter_outputFile, "# Tempo residuo: %12g\n", tempoResiduo);
    fprintf(parameter_outputFile, "#\n");
    fprintf(parameter_outputFile, "n: %d\n", parameter_n);
    fprintf(parameter_outputFile, "x:\n");
    fprnVetor(parameter_outputFile ,x, parameter_n);

    liberaSisDiag(slOtmSimetrico);
    liberaSisDiag(slOtm);
    free(x);
    free(maxNormaInter);

    free(MARKER_NAME_op1);
    free(MARKER_NAME_op2);
    LIKWID_MARKER_CLOSE;

    return status;
}