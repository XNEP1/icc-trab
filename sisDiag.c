#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sisDiag.h"

// Retorna a quantidade de numeros não numeros que existem na linha 'i'.
int qntNumLinha(int i, int n, int k){
    int q = (k-1)/2;

    if(i < q)
        return k - (q - i);
    else if(i > n-q-1)
        return k - (q - (n-i-1));
    else
        return k;
}

// Retorna a quantidade de números não nulos numa matriz k-diagonal
int qntNumMatriz(int n, int k){
    return n + (((4*n - k - 1)*(k-1)) / 4) + qntNumLinha(n - 1, n, k);
}

// Alocaçao e desalocação de matrizes
SisDiag_t* alocaSisDiag (unsigned int n, unsigned int k){
    SisDiag_t *SD = (SisDiag_t*) malloc(sizeof(SisDiag_t));

    if ( SD == NULL ) {
        return NULL;
    }

    SD->n = n;
    SD->k = k;

    SD->A = (real_t **) malloc(n * sizeof(real_t *));
    SD->b = (real_t *) calloc(n, sizeof(real_t));

    if (!(SD->A) || !(SD->b)) {
        liberaSisDiag(SD);
        return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SD->A[0] = (real_t *) calloc(qntNumMatriz(n, k) , sizeof(real_t));
    if (!(SD->A[0])) {
        liberaSisDiag(SD);
        return NULL;
    }

    int q = (k-1)/2; // quantidade de diagonais acima(ou abaixo) da diagonal principal.

    for (int i=1; i <= q; ++i) {
        SD->A[i] = SD->A[i-1]+ qntNumLinha(i-1, n, k);
    }
    
    SD->A[q+1] = SD->A[q] + (qntNumLinha(q, n, k) - 1);

    for (int i=q+2; i < n; ++i) {
        SD->A[i] = SD->A[i-1] + (i-q-1) + qntNumLinha(i-1, n, k) - (i-q);
    }

    return (SD);
}

void liberaSisDiag (SisDiag_t *SD){
  if (SD) {
    if (SD->A) {
      if (SD->A[0]) free (SD->A[0]);
    free(SD->A);
    }
    
    if (SD->b) free(SD->b);

    free(SD);
  }
}

// Leitura e impressão de sistemas lineares
void prnSisDiag (SisDiag_t *SD){

    unsigned int n = SD->n;

  for(int i=0; i < n; ++i) {
    printf("\n  ");
    for(int j=0; j < n; ++j){
        if(IS_ZERO(SD,i,j))
            printf("%10g", 0.0);
        else
            printf ("%10g", SD->A[i][j]);
    }

    printf ("   |   %g", SD->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  int i;

    for(i=0; i < n; ++i)
        printf ("%.15g ", v[i]);
    
    printf("\n\n");

}

void fprnVetor (FILE *file,real_t *v, unsigned int n){
  int i;

    for(i=0; i < n; ++i)
        fprintf (file ,"%.15g ", v[i]);
    
    fprintf(file, "\n\n");
}

/***********************
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
 * k: numero de diagonais da matriz A
 ***********************/
static inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}

/***********************
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * k: numero de diagonais da matriz A
 ***********************/
static inline double generateRandomB( unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k<<2) * (double)rand() * invRandMax;
}

void gerarCoeficientesSD(SisDiag_t *SD){
    int i, j;
    unsigned int n = SD->n;

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){

            if(IS_ZERO(SD, i, j))
                continue; // Fora da banda de diagonais.

            SD->A[i][j] = generateRandomA(i, j, SD->k);

        }
        
        SD->b[i] = generateRandomB(SD->k);
    }
}