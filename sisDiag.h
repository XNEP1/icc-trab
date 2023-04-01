#ifndef __SISDIAG_H__
#define __SISDIAG_H__

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

/* Copia os valores de source em dest. */
/* dest já deve estar alocado e deve ter o mesmo n de source. */
int copiaSisDiag(SisDiag_t *source, SisDiag_t *dest);

// Leitura e impressão de sistemas lineares
void prnSisDiag (SisDiag_t *SD);

void prnVetor (real_t *v, unsigned int n);
void fprnVetor (FILE *file,real_t *v, unsigned int n);

void gerarCoeficientesSD(SisDiag_t *SD);

#endif // __SISDIAG_H__