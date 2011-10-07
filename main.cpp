#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <cassert>

#define MAXM 100
#define MAXV 100

// #define DEBUG 1
// #define DEBUGLAMBDA 1
// #define DEBUGSUBGRADIENT 1
#define EPSILON 10e-9

using namespace std;

int m, k;
int bh[MAXM], bv[MAXM]; // capacidade de memória do satélite para as faixas
// int r[MAXM][MAXM]; // ganho pela captura do fragmento i j
int **r;
// int mh[MAXM][MAXM], mv[MAXM][MAXM]; // memória necessária para capturar fragmento ij em cada observacao.
int **mh, **mv;
int v; // numero de vertices do poligono
int vx[MAXV], vy[MAXV];
// double lambda[MAXM][MAXM];
double **lambda;
int xh[MAXM][MAXM], xv[MAXM][MAXM];
double **dp;
int maxiter;

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

    // aloca matrizes para guardar o ganho e pesos de cada shard.
    r = new int*[m+1];
    mh = new int*[m+1];
    mv = new int*[m+1];
    for (int i = 0; i <= m; i++) {
      r[i] = new int[m+1];
      mh[i] = new int[m+1];
      mv[i] = new int[m+1];
    }

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
double solveLRi(int f, int W) {
  // aloca matriz de programacao dinamica da mochila
  dp = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    dp[i] = new double[W+1]; // matriz m X W (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i <= m; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i <= W; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= W; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mh[f][i] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mh[f][i])] + (r[f][i] - lambda[f][i]));
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  double res = dp[m][W];

  int i = m, j = W;
  while (i) {
    if (dp[i][j] == dp[i-1][j]) {
      xh[f][i] = 0;
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mh[f][i])] + (r[f][i] - lambda[f][i]))) {
      xh[f][i] = 1;
      j = max(0, j - mh[f][i]);
      i--;
    }
  }

  #ifdef DEBUG
  printSolution(f, true);
  #endif

  // desaloca matriz da mochila
  for (int i = 0; i <= m; i++) {
    delete[] dp[i];
  }
  delete[] dp;

  return res;
}

// Resolve a mochila da faixa f
// As faixas f são verticais
double solveLRj(int f, int W) {
  // aloca matriz de programacao dinamica da mochila
  dp = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    dp[i] = new double[W+1]; // matriz m X W (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i <= m; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i <= W; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= W; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mv[i][f] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mv[i][f])] + (r[i][f] - lambda[i][f]));
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  double res = dp[m][W];

  int i = m, j = W;
  while (i) {
    if (dp[i][j] == dp[i-1][j]) {
      xv[i][f] = 0;
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mv[i][f])] + (r[i][f] - lambda[i][f]))) {
      xv[i][f] = 1;
      j = max(0, j - mv[i][f]);
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

  return res;
}

double solveLR() {
  double z = 0;
  for (int i = 1; i <= m; i++) {
    // printf("%d\r", i);
    // fflush(stdout);
    z += solveLRi(i, bh[i]);
    z += solveLRj(i, bv[i]);
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

// Resolve a mochila da faixa f considerando apenas os itens disponiveis após fixação
// As faixas f são horizontais
double mochilaHorizontal(int f, int W, vector <int> disponiveis) {
  int sz = disponiveis.size();

  if (W <= 0) {
    for (int i = 1; i <= sz; i++) {
      xh[f][disponiveis[i-1]] = 0;
    }
    return 0;
  }

  // aloca matriz de programacao dinamica da mochila
  dp = new double*[sz+1];
  for (int i = 0; i <= sz; i++) {
    dp[i] = new double[W+1]; // matriz sz X W (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i <= sz; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i <= W; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= sz; i++) {
    for (int j = 1; j <= W; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mh[f][disponiveis[i-1]] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mh[f][disponiveis[i-1]])] + r[f][disponiveis[i-1]]);
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }


  double res = dp[sz][W];

  int i = sz, j = W;
  while (i) {
    if (dp[i][j] == dp[i-1][j]) {
      xh[f][disponiveis[i-1]] = 0;
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mh[f][disponiveis[i-1]])] + r[f][disponiveis[i-1]])) {
      xh[f][disponiveis[i-1]] = 1;
      j = max(0, j - mh[f][disponiveis[i-1]]);
      i--;
    }
  }

  #ifdef DEBUG
  printSolution(f, true);
  #endif

  // desaloca matriz da mochila
  for (int i = 0; i <= sz; i++) {
    delete[] dp[i];
  }
  delete[] dp;

  return res;
}

// Resolve a mochila da faixa f considerando apenas os itens disponiveis após fixação
// As faixas f são horizontais
double mochilaVertical(int f, int W, vector <int> disponiveis) {
  int sz = disponiveis.size();

  if (W <= 0) {
    for (int i = 1; i <= sz; i++) {
      xv[disponiveis[i-1]][f] = 0;
    }
    return 0;
  }

  // aloca matriz de programacao dinamica da mochila
  dp = new double*[sz+1];
  for (int i = 0; i <= sz; i++) {
    dp[i] = new double[W+1]; // matriz sz X W (capacidade maxima da mochila)
  }
  
  // zerando primeira coluna
  for (int i = 0; i <= sz; i++) { 
    dp[i][0] = 0;
  }
  // zerando primeira linha
  for (int i = 0; i <= W; i++) {
    dp[0][i] = 0;
  }

  for (int i = 1; i <= sz; i++) {
    for (int j = 1; j <= W; j++) {
      // dp[i][j] = máximo entre {levar o item atual + maximo com a capacidade restante} e {sem o item atual}
      if (mv[disponiveis[i-1]][f] <= j) {
        dp[i][j] = max(dp[i-1][j], dp[i-1][max(0, j - mv[disponiveis[i-1]][f])] + r[disponiveis[i-1]][f]);
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  int res = dp[sz][W];

  int i = sz, j = W;
  while (i) {
    if (dp[i][j] == dp[i-1][j]) {
      xv[disponiveis[i-1]][f] = 0;
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mv[disponiveis[i-1]][f])] + r[disponiveis[i-1]][f])) {
      xv[disponiveis[i-1]][f] = 1;
      j = max(0, j - mv[disponiveis[i-1]][f]);
      i--;
    }
  }

  #ifdef DEBUG
  printSolution(f, true);
  #endif

  // desaloca matriz da mochila
  for (int i = 0; i <= sz; i++) {
    delete[] dp[i];
  }
  delete[] dp;

  return res;
}


double heuristicLB() {
  int lb = 0;
  vector <int> disponiveis;
  int valorFixo;
  int pesoRestante;

  for (int i = 1; i <= m; i++) {
    disponiveis.clear();
    pesoRestante = bh[i];
    valorFixo = 0;
    for (int j = 1; j <= m; j++) {
      if (xh[i][j] + xv[i][j] != 1) {
        disponiveis.push_back(j);
        pesoRestante -= mh[i][j];
      }
      else if (xh[i][j] == 1){
        valorFixo += r[i][j];
      }
    }
    lb += mochilaHorizontal(i, pesoRestante, disponiveis) + valorFixo;
  }
  // verificar quais os novos selecionados e fixar.
  for (int j = 1; j <= m; j++) {
    disponiveis.clear();
    pesoRestante = bv[j];
    valorFixo = 0;
    for (int i = 1; i <= m; i++) {
      if (xh[i][j] == 0 && xv[i][j] == 0) { // só estão disponíveis os que ainda não foram selecionadas para a passada vertical
        disponiveis.push_back(i);
        pesoRestante -= mv[i][j];
      }
      else if (xh[i][j] == 0 && xv[i][j] == 1) {
        valorFixo += r[i][j];
      }
      else if (xh[i][j] == 1 && xv[i][j] == 1) { // se já foi selecionado pelo horizontal
        xv[i][j] = 0; // então sai desta mochila
      }
    }
    lb += mochilaVertical(j, pesoRestante, disponiveis) + valorFixo;
  }

  int totalGain = 0;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      if (xh[i][j] == 1)
        totalGain += r[i][j];
      if (xv[i][j] == 1) 
        totalGain += r[i][j];
      if (xh[i][j] == 1 && xv[i][j] == 1) {
        cout << "FOOOOODEUUUUUU" << endl;
        cout << "FOOOOODEUUUUUU" << endl;
        cout << "FOOOOODEUUUUUU" << endl;
      }
    }
  }
  
  // cout << " ----------------------------> " << totalGain << " " << lb << endl;

  return lb;
}

void subgradientOpt(double pi) {
  
  double T;
  double zLB, zLB_current; // uma solução viável via heurística.
  double zUP = 123123123; // limitante dual (superior)
  double zUP_current; // limitante dual (superior) corrente
  int G[MAXM][MAXM];
  int sumSquaredG;

  // aloca matriz com os lambdas
  lambda = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    lambda[i] = new double[m+1]; 
  }

  initializeLambdas(); 

  zLB = 0; // ainda não temos um bom lower bound, então fica como 0
  while (maxiter--) {
    zUP_current = solveLR(); // resolve relaxação lagrangeana
    // cout << zUP_current << endl;
    if (zUP_current < zUP) {
      zUP = zUP_current;
      cout << "Novo limitante dual: " << setiosflags(ios::fixed) << setprecision(6) << zUP << endl;
    }
    updateG(G, &sumSquaredG); // atualiza os subgradientes e calcula a soma dos Gi^2
    T = pi*(zUP - 0) / sumSquaredG; // atualiza o tamanho do passo
    updateLambdas(T, G);
    pi = 0.99 * pi;
    zLB_current = heuristicLB();
    if (zLB_current > zLB) {
      zLB = zLB_current;
      cout << "Novo limitante primal: " << zLB << endl;
    }
  }

  // desaloca matriz de lambdas
  for (int i = 0; i <= m; i++) {
    delete[] lambda[i];
  }
  delete[] lambda;
}

int main(int argc, char **argv) {

  if (argc < 3) {
    cout << "USAGE: " << argv[0] << " ENTRADA SAIDA" << endl;
    return 1;
  }

  // tenta ler entrada
  bool feedback;
  feedback = readInput(argv[1]);
  if (!feedback) {
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
    param >> maxiter;
  }
  catch (ios_base::failure) {
    cout << "Erro na leitura do arquivo de parâmetros." << endl;
    return false;
  }


  // double z = solveLR();
  // cout << endl << z << endl;

  subgradientOpt(pi);

  // desaloca matrizes que guardam ganhos e pesos
  for (int i = 0; i <= m; i++) {
    delete[] r[i];
    delete[] mv[i];
    delete[] mh[i];
  }
  delete[] r;
  delete[] mv;
  delete[] mh;


  return 0;
}
