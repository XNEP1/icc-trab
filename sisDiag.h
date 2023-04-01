#ifndef __SISDIAG_H__
#define __SISDIAG_H__

// retorna 1 se a for menor que b. 0 caso contrario.
#define MENOR(a, b) (unsigned int)(((unsigned)((a) - (b)) >> (sizeof(int) * CHAR_BIT - 1)) & 1)

// 1 se nesses indeces há zero ou 0 se nesses indices não é de zero.
#define IS_ZERO(SL, i, j) (ABS((j)-(i)) > (SL->k-1)/2)

#define GET_A(SL, i, j) ((IS_ZERO(SL, i, j)) ? 0.0 : SL->A[i][j])

// bit mais significativo
#define MSB (sizeof(int) * CHAR_BIT - 1)

// 1 se é positivo. 0 se é negativo
#define IS_POSITIVE(a) (((unsigned int)( (a) >> MSB) & 1) ^ 1)

// Estrutura para definiçao de um sistema linear com coeficientes tridiagonais (p = q = 1)
typedef struct {
  real_t **A; // coeficientes
  real_t *b; // termos independentes
  unsigned int n; // tamanho do SL
  unsigned int k; // banda do SL k-diagonal
} SisDiag_t;

// Alocaçao e desalocação de matrizes
SisDiag_t* alocaSisDiag (unsigned int n, unsigned int k);
void liberaSisDiag (SisDiag_t *SD);

// Leitura e impressão de sistemas lineares
void prnSisDiag (SisDiag_t *SD);

void prnVetor (real_t *v, unsigned int n);
void fprnVetor (FILE *file,real_t *v, unsigned int n);

void gerarCoeficientesSD(SisDiag_t *SD);

#endif // __SISDIAG_H__