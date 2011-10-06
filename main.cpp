#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cmath>

#define MAXM 1000
#define MAXV 1000

#define DEBUG 1
#define DEBUGLAMBDA 1
#define DEBUGSUBGRADIENT 1
#define EPSILON 10e-9

using namespace std;

int m, k;
int bh[MAXM], bv[MAXM]; // capacidade de memória do satélite para as faixas
int r[MAXM][MAXM]; // ganho pela captura do fragmento i j
int mh[MAXM][MAXM], mv[MAXM][MAXM]; // memória necessária para capturar fragmento ij em cada observacao.
int v; // numero de vertices do poligono
int vx[MAXV], vy[MAXV];
double lambda[MAXM][MAXM];
int xh[MAXM][MAXM], xv[MAXM][MAXM];
double res;
double **dp;


bool equal(double a, double b) {
  return fabs(a-b) < EPSILON;
}

void initializeLambdas() {
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= m; j++) {
      lambda[i][j] = 0;
    }
  }
}

bool readInput(char *f) {
  ifstream input(f);
  
  if (!input.is_open())
    return false;

  try {
    input >> m; // número de faixas
    int i, bi, j;
    for (int p = 0; p < m; p++) { // lê faixas horizontais e seus bi;
      input >> i;
      input >> bi;
      bh[i] = bi;
    } 
    for (int p = 0; p < m; p++) { // lê faixas verticais e seus bi;
      input >> i;
      input >> bi;
      bv[i] = bi;
    }
    input >> k;
    for (int p = 0; p < k; p++) {
      input >> i >> j;
      input >> r[i][j];
      input >> mh[i][j];
      input >> mv[i][j];
    }
    input >> v;
    for (int p = 0; p < v; p++) {
      input >> vx[p];
      input >> vy[p];
    }
    input.close(); 
  }
  catch (ios_base::failure) {
    cout << "Erro na leitura do arquivo depois de aberto." << endl;
    return false;
  }
  return true;
}

void printInput() {
  cout << m << endl;
  for (int p = 0; p < m; p++) { 
    cout << p+1<< " " << bh[p+1] << " ";
  }
  cout << endl;
  for (int p = 0; p < m; p++) { 
    cout << p+1 << " " << bv[p+1] << " ";
  }
  cout << endl << k << endl;
}


void printSolution(int f, bool horiz) {

  cout << endl << "Faixa " << f << " ";
  if (horiz) {
    cout << "H" << endl;
    // imprime matriz da dp da mochila
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= bh[f]; j++) {
        cout << setw(2) << dp[i][j] << " ";
      }
      cout << endl;
    }

    // imprime tabela com os items selecionados
    cout << endl << endl << "   ";
    for (int j = 1; j <= m; j++) {
      cout <<  setw(2) << j << " ";
    }
    cout << endl;

    for (int i = 1; i <= m; i++) {
      cout << setw(2) << i << " ";
      for (int j = 1; j <= m; j++) {
        cout << setw(2) << xh[i][j] << " ";
      }
      cout << endl;
    }
  }
  else {
    cout << "V" << endl;
    // imprime matriz da dp da mochila
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= bv[f]; j++) {
        cout << setw(2) << dp[i][j] << " ";
      }
      cout << endl;
    }

    // imprime tabela com os items selecionados
    cout << endl << endl << "   ";
    for (int j = 1; j <= m; j++) {
      cout <<  setw(2) << j << " ";
    }
    cout << endl;

    for (int i = 1; i <= m; i++) {
      cout << setw(2) << i << " ";
      for (int j = 1; j <= m; j++) {
        cout << setw(2) << xv[i][j] << " ";
      }
      cout << endl;
    }

  }

}

// Resolve a mochila da faixa f
// As faixas f são horizontais
void solveLRi(int f) {
  // aloca matriz de programacao dinamica da mochila
  dp = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    dp[i] = new double[bh[f]+1]; // matriz m X bh[f] (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i <= m; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i <= bh[f]; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= bh[f]; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mh[f][i-1] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mh[f][i-1])] + (r[f][i-1] - lambda[f][i]));
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  res = dp[m][bh[f]];

  int i = m, j = bh[f];
  while (i) {
    if (dp[i][j] == dp[i-1][j]) {
      xh[f][i-1] = 0;
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mh[f][i-1])] + (r[f][i-1] - lambda[f][i]))) {
      xh[f][i-1] = 1;
      j = max(0, j - mh[f][i-1]);
      i--;
    }
    // else {
    //   xh[f][i-1] = 1;
    //   j = max(0, j - mh[f][i-1]);
    //   i--;
    // }
  }

  #ifdef DEBUG
  printSolution(f, true);
  #endif

  // desaloca matriz da mochila
  for (int i = 0; i <= m; i++) {
    delete[] dp[i];
  }
  delete[] dp;
}

// Resolve a mochila da faixa f
// As faixas f são verticais
void solveLRj(int f) {
  // aloca matriz de programacao dinamica da mochila
  dp = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    dp[i] = new double[bv[f]+1]; // matriz m X bv[f] (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i <= m; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i <= bv[f]; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= bv[f]; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mv[i-1][f] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mv[i-1][f])] + (r[i-1][f] - lambda[i][f]));
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  res = dp[m][bv[f]];

  int i = m, j = bv[f];
  while (i) {
    if (dp[i][j] == dp[i-1][j]) {
      xv[i-1][f] = 0;
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mv[i-1][f])] + (r[i-1][f] - lambda[i][f]))) {
      xv[i-1][f] = 1;
      j = max(0, j - mv[i-1][f]);
      i--;
    }
  }

  #ifdef DEBUG
  printSolution(f, false);
  #endif

  // desaloca matriz da mochila
  for (int i = 0; i <= m; i++) {
    delete[] dp[i];
  }
  delete[] dp;
}

double solveLR() {
  double z;
  for (int i = 1; i <= m; i++) {
    // printf("%d\r", i);
    // fflush(stdout);
    solveLRi(i);
    z += res;
    solveLRj(i);
    z += res;
  }
  // cout << endl;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      z += lambda[i][j];
    }
  }
  return z;
}

void updateG(int G[MAXM][MAXM], int *sumSquaredG) {
  *sumSquaredG = 0;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      G[i][j] = 1 - xh[i][j] - xv[i][j];
      *sumSquaredG += G[i][j] * G[i][j];
      #ifdef DEBUGSUBGRADIENT
      cout << setw(3) << G[i][j] << " ";
      #endif
    }
    #ifdef DEBUGSUBGRADIENT
    cout << endl;
    #endif
  }
}

void updateLambdas(double T, int G[MAXM][MAXM]) {
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      lambda[i][j] = max(0.0, lambda[i][j] - T*G[i][j]);
      #ifdef DEBUGLAMBDA
      cout << setw(3) << lambda[i][j] << " ";
      #endif
    }
    #ifdef DEBUGLAMBDA
    cout << endl;
    #endif
  }

}

void subgradientOpt(double pi) {
  
  double T;
  double zLB; // uma solução viável via heurística.
  double zUP = 123123123; // limitante dual (superior)
  double zUP_current; // limitante dual (superior) corrente
  int G[MAXM][MAXM];
  int sumSquaredG;

  zLB = 0;
  int iter = 200;
  while (iter--) {
    zUP_current = solveLR(); // resolve relaxação lagrangeana
    cout << zUP_current << endl;
    if (zUP_current < zUP) {
      zUP = zUP_current;
      cout << "Novo limitante dual: " << zUP << endl;
    }
    updateG(G, &sumSquaredG); // atualiza os subgradientes e calcula a soma dos Gi^2
    T = pi*(zUP - zLB) / sumSquaredG; // atualiza o tamanho do passo
    updateLambdas(T, G);
  }

}

int main(int argc, char **argv) {

  if (argc < 3) {
    cout << "USAGE: " << argv[0] << " ENTRADA SAIDA" << endl;
    return 1;
  }

  // tenta ler entrada
  bool r;
  r = readInput(argv[1]);
  if (!r) {
    cout << "Erro ao ler entrada!" << endl;
    return 1;
  }

  double pi;
  ifstream param("param");
  if (!param.is_open()) {
    cout << "Não foi possível abrir o arquivo de parâmetros." << endl;
    return false;
  }
  try{
    param >> pi;
  }
  catch (ios_base::failure) {
    cout << "Erro na leitura do arquivo de parâmetros." << endl;
    return false;
  }

  initializeLambdas();

  double z = solveLR();
  cout << endl << z << endl;

  // subgradientOpt(pi);

  return 0;
}
