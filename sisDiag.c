#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sisDiag.h"

// Alocaçao e desalocação de matrizes
SisDiag_t* alocaSisDiag (unsigned int n, unsigned int k){
    SisDiag_t *SD = (SisDiag_t*) malloc(sizeof(SisDiag_t));

    if ( SD ) {

    SD->n = n;
    SD->k = k;

    SD->A = (real_t **) malloc(n * sizeof(real_t *));
    SD->b = (real_t *) calloc(n, sizeof(real_t));

    if (!(SD->A) || !(SD->b)) {
        liberaSisDiag(SD);
        return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SD->A[0] = (real_t *) calloc(n * k , sizeof(real_t));
    if (!(SD->A[0])) {
        liberaSisDiag(SD);
        return NULL;
    }

    for (int i=1; i < n; ++i) {
        SD->A[i] = SD->A[i-1]+k;
    }
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

/* Copia os valores de source em dest. */
/* dest já deve estar alocado e deve ter o mesmo n de source. */
int copiaSisDiag(SisDiag_t *source, SisDiag_t *dest){
    int i, j;
    unsigned int n = source->n;
    unsigned int k = source->k;
    if(k == 0)
    k = n;

    if(source->n != dest->n)
    return 0;

    for(i=0; i<n; i++){
    for(j=0; j<k; j++)
        dest->A[i][j] = source->A[i][j];

    dest->b[i] = source->b[i];
    }

    return 1;
}

// Leitura e impressão de sistemas lineares
void prnSisDiag (SisDiag_t *SD){

    unsigned int n = SD->n;
    unsigned int k = SD->k;

  for(int i=0; i < n; ++i) {
    printf("\n  ");
    for(int j=0; j < k; ++j)
      printf ("%10g", SD->A[i][j]);
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
    int q = SD->k / 2;

    for(i=0; i<n; i++)
        for(j=0; j<n; j++){
            if( (j>=(i-q)) && (j<=(i+q))){
                SD->A[i][j - i + q] = generateRandomA(i, j, SD->k);
            }

            SD->b[i] = generateRandomB(SD->k);
        }
}