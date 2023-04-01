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

int main(int argc, char **argv){

    srand(20222);
    LIKWID_MARKER_INIT;

    int parameter_n;
    int parameter_k;
    int parameter_p;
    int parameter_i;
    real_t parameter_e = -1.0;
    FILE *parameter_outputFile; 

    SisDiag_t *sl;
    real_t *x; // resultado
    real_t *maxNormaInter; // inter k
    real_t normaResiduo;
    rtime_t tempoPC = 0.0;
    rtime_t tempoInter = 0.0;
    rtime_t tempoResiduo = 0.0;
    rtime_t tempoTotal = timestamp();
    rtime_t tempoAux;
    int i, j;
    int interacoes;
    int status = 0;
    int precondicionador;

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

    if( parameter_k > (2*parameter_n - 1)){
        fprintf(stderr, "O valor de K está muito alto para o esse valor de n. (N=%d)(K=%d)(Max=%d)", parameter_n , parameter_k, (2*parameter_n - 1) );
        exit(-1);
    }

    sl = alocaSisDiag(parameter_n, parameter_k);
    x = (real_t *) calloc(parameter_n , sizeof(real_t)); // Todos os valores de x começam sendo 0.0
    maxNormaInter = (real_t *) calloc(parameter_i, sizeof(real_t));

    if(!(sl) || !(x) || !(maxNormaInter)){
        perror("Erro de alocação de memoria.\n");
        exit(-2);
    }

    gerarCoeficientesSD(sl);

    if(parameter_p == 0)
        precondicionador = PRECONDIONADOR_NENHUM;
    else
        precondicionador = PRECONDIONADOR_JACOBI;

    string_t MARKER_NAME_op1 = markerName("conjGradient", parameter_n);
    LIKWID_MARKER_START(MARKER_NAME_op1);
    interacoes = conjGradient(sl, x, parameter_e, parameter_i, precondicionador ,&tempoInter, &tempoPC,maxNormaInter);
    LIKWID_MARKER_STOP(MARKER_NAME_op1);

    string_t MARKER_NAME_op2 = markerName("calResiduo", parameter_n);
    real_t *r = malloc(parameter_n * sizeof(real_t));
    LIKWID_MARKER_START(MARKER_NAME_op2);
    calResiduo(sl, x, r);
    LIKWID_MARKER_STOP(MARKER_NAME_op2);

    tempoAux  = timestamp();
    normaResiduo = normaL2Residuo(sl, x);
    tempoResiduo = timestamp() - tempoAux;

    if(normaResiduo >= 1.0){
        fprintf(stderr, "n: %d > Residuo muito alto ( %10g ): o metodo pode não ter encontrado uma solução.\n", parameter_n, normaResiduo);
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
    fprintf(parameter_outputFile, "# Tempo Total: %12g\n", timestamp() - tempoTotal);
    fprintf(parameter_outputFile, "#\n");
    fprintf(parameter_outputFile, "n: %d\n", parameter_n);
    fprintf(parameter_outputFile, "x:\n");
    fprnVetor(parameter_outputFile ,x, parameter_n);

    liberaSisDiag(sl);
    free(x);
    free(maxNormaInter);

    free(MARKER_NAME_op1);
    free(MARKER_NAME_op2);
    LIKWID_MARKER_CLOSE;

    return status;
}