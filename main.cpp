#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>

#define MAXM 10000
#define MAXV 10000

using namespace std;

int m, k;
int bh[MAXM], bv[MAXM]; // capacidade de memória do satélite para as faixas
int r[MAXM][MAXM]; // ganho pela captura do fragmento i j
int mh[MAXM][MAXM], mv[MAXM][MAXM]; // memória necessária para capturar fragmento ij em cada observacao.
int v; // numero de vertices do poligono
int vx[MAXV], vy[MAXV];
double lambda[MAXM][MAXM];
double res;

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

// Resolve a mochila da faixa f
// As faixas f são horizontais
void solveLRi(int f) {
  // aloca matriz de programacao dinamica da mochila
  double **dp = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    dp[i] = new double[bh[f]+1]; // matriz m X bh[f] (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i < m; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i < bh[f]; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= bh[f]; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mh[f][i-1] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mh[f][i-1])] + r[f][i-1]);
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  res = dp[m][bh[f]];

  // imprime matriz da dp da mochila
  // for (int i = 0; i <= m; i++) {
  //   for (int j = 0; j <= bh[f]; j++) {
  //     cout << setw(2) << dp[i][j] << " ";
  //   }
  //   cout << endl;
  // }

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
  double **dp = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    dp[i] = new double[bv[f]+1]; // matriz m X bv[f] (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i < m; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i < bv[f]; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= bv[f]; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mv[i-1][f] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mv[i-1][f])] + r[i-1][f]);
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  res = dp[m][bv[f]];

  // imprime matriz da dp da mochila
  // for (int i = 0; i <= m; i++) {
  //   for (int j = 0; j <= bv[f]; j++) {
  //     cout << setw(2) << dp[i][j] << " ";
  //   }
  //   cout << endl;
  // }

  // desaloca matriz da mochila
  for (int i = 0; i <= m; i++) {
    delete[] dp[i];
  }
  delete[] dp;
}

double solveLR() {
  double z;
  for (int i = 1; i <= m; i++) {
    printf("%d\r", i);
    fflush(stdout);
    solveLRi(i);
    z += res;
    solveLRj(i);
    z += res;
  }
  cout << endl;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      z += lambda[i][j];
    }
  }
  return z;
}

void subgradientOpt(double pi) {
  
  double T;
  double zLB; // uma solução viável via heurística.
  double zUP; // limitante dual (superior)
  int Gh[MAXM], Gv[MAXM];
  int sumSquaredG;

  int iter = 100;
  while (iter--) {
    solveLR(); // resolve relaxação lagrangeana
    updateG(Gh, Gv, &sumSquaredG); // atualiza os subgradientes e calcula a soma dos Gi^2
    T = pi*(zUP - zLB) / sumSquaredG; // atualiza o tamanho do passo
    updateLambdas(T, Gh, Gv);
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

  double z = solveLR();
  cout << z << endl;

  return 0;
}
