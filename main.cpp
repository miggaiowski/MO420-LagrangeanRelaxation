#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include <signal.h>
#include <cfloat>

#define NDEBUG 
#include <cassert>


#define MAXM 1000
#define MAXV 1000
#define NUM_THREADS 16

// #define DEBUG 1
// #define DEBUGLAMBDA 1
// #define DEBUGSUBGRADIENT 1
#define EPSILON 10e-7

using namespace std;

int m, k;
int bh[MAXM], bv[MAXM]; // capacidade de memória do satélite para as faixas
int **r; // ganho pela captura do fragmento i j
int **mh, **mv; // memória necessária para capturar fragmento ij em cada observacao.
int v; // numero de vertices do poligono
int vx[MAXV], vy[MAXV];
double **lambda;
char xh[MAXM][MAXM], xv[MAXM][MAXM];
char temp_xh[MAXM][MAXM], temp_xv[MAXM][MAXM];
char temp_xh2[MAXM][MAXM], temp_xv2[MAXM][MAXM];
char best_xh[MAXM][MAXM], best_xv[MAXM][MAXM];
int maxiter;
double zUP, zLB;
int nSemMelhora = 50;
time_t startTime, currentTime, maxTime;
double menorPi;
int iter = 0;
time_t primalTime, dualTime;

void quit(int pa) {
  cerr << "Recebi SIGINT, saindo depois desta iteração." << endl;
  cerr << "Iterações executadas: " << iter << endl;
  iter = maxiter;
}

bool equal(double a, double b) {
  return fabs(a-b) < EPSILON;
}

void copiaMatrizesX() {
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= m; j++) {
      best_xh[i][j] = xh[i][j];
      best_xv[i][j] = xv[i][j];
    }
  }
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

  // cout << endl << "Faixa " << f << " ";
  // if (horiz) {
  //   cout << "H" << endl;
  //   // imprime matriz da dp da mochila
  //   for (int i = 0; i <= m; i++) {
  //     for (int j = 0; j <= bh[f]; j++) {
  //       cout << setw(2) << dp[i][j] << " ";
  //     }
  //     cout << endl;
  //   }

  //   // imprime tabela com os items selecionados
  //   cout << endl << endl << "   ";
  //   for (int j = 1; j <= m; j++) {
  //     cout <<  setw(2) << j << " ";
  //   }
  //   cout << endl;

  //   for (int i = 1; i <= m; i++) {
  //     cout << setw(2) << i << " ";
  //     for (int j = 1; j <= m; j++) {
  //       cout << setw(2) << xh[i][j] << " ";
  //     }
  //     cout << endl;
  //   }
  // }
  // else {
  //   cout << "V" << endl;
  //   // imprime matriz da dp da mochila
  //   for (int i = 0; i <= m; i++) {
  //     for (int j = 0; j <= bv[f]; j++) {
  //       cout << setw(2) << dp[i][j] << " ";
  //     }
  //     cout << endl;
  //   }

  //   // imprime tabela com os items selecionados
  //   cout << endl << endl << "   ";
  //   for (int j = 1; j <= m; j++) {
  //     cout <<  setw(2) << j << " ";
  //   }
  //   cout << endl;

  //   for (int i = 1; i <= m; i++) {
  //     cout << setw(2) << i << " ";
  //     for (int j = 1; j <= m; j++) {
  //       cout << setw(2) << xv[i][j] << " ";
  //     }
  //     cout << endl;
  //   }

  // }

}

// Resolve a mochila da faixa f
// As faixas f são horizontais
double solveLRi(int f, int W) {
  double **dp;

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
        dp[i][j] = max(dp[i-1][j], dp[i-1][j - mh[f][i]] + (r[f][i] - lambda[f][i]));
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
      #pragma omp critical
      {
      xh[f][i] = 0;
      }
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mh[f][i])] + (r[f][i] - lambda[f][i]))) {
      #pragma omp critical
      {
      xh[f][i] = 1;
      }
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
  double **dp;

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
        dp[i][j] = max(dp[i-1][j], dp[i-1][j - mv[i][f]] + (r[i][f] - lambda[i][f]));
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
      #pragma omp critical
      {
      xv[i][f] = 0;
      }
      i--;
    }
    else if (equal(dp[i][j], dp[i-1][max(0, j - mv[i][f])] + (r[i][f] - lambda[i][f]))) {
      #pragma omp critical
      {
      xv[i][f] = 1;
      }
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
  #pragma omp parallel for reduction(+:z)
  for (int i = 1; i <= m; i++) {
    // double z1, z2;
    // printf("%d\r", i);
    // fflush(stdout); 
    z += solveLRi(i, bh[i]);
    z += solveLRj(i, bv[i]);
    // #pragma omp atomic
    // z += z1+z2;
  }
  
  // cout << endl;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      z += lambda[i][j];
    }
  }
  return z;
}


// Resolve a mochila da faixa f considerando apenas os itens disponiveis após fixação
// As faixas f são horizontais
double mochilaHorizontal(int f, int W, vector <int> disponiveis) {
  double **dp;

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
        dp[i][j] = max(dp[i-1][j], dp[i-1][j - mh[f][disponiveis[i-1]]] + r[f][disponiveis[i-1]]);
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }


  double res = dp[sz][W];

  int i = sz, j = W;
  while (i) {
    if (equal(dp[i][j], dp[i-1][j])) {
      xh[f][disponiveis[i-1]] = 0;
      i--;
    }
    else if (equal(dp[i][j], 
                   dp[i-1][max(0, j - mh[f][disponiveis[i-1]])] + r[f][disponiveis[i-1]])) {
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
  double **dp;

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
        dp[i][j] = max(dp[i-1][j], dp[i-1][j - mv[disponiveis[i-1]][f]] + r[disponiveis[i-1]][f]);
      }
      else {
        dp[i][j] = dp[i-1][j];
      }
    }
  }

  int res = dp[sz][W];

  int i = sz, j = W;
  while (i) {
    if (equal(dp[i][j], dp[i-1][j])) {
      xv[disponiveis[i-1]][f] = 0;
      i--;
    }
    else if (equal(dp[i][j], 
                   dp[i-1][max(0, j - mv[disponiveis[i-1]][f])] + r[disponiveis[i-1]][f])) {
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

void copyTempX(int qual) {
  if (qual == 1) {
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= m; j++) {
        temp_xv[i][j] = xv[i][j];
        temp_xh[i][j] = xh[i][j];
      }
    }
  }
  else {
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= m; j++) {
        temp_xv2[i][j] = xv[i][j];
        temp_xh2[i][j] = xh[i][j];
      }
    }
  }
}

void copyTempXBack(int qual) {
  if (qual == 1) {
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= m; j++) {
        xv[i][j] = temp_xv[i][j];
        xh[i][j] = temp_xh[i][j];
      }
    }
  }
  else {
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= m; j++) {
        xv[i][j] = temp_xv2[i][j];
        xh[i][j] = temp_xh2[i][j];
      }
    }
  }
}

double heuristicLB() {
  int lb = 0, lbh = 0, lbv = 0;
  int lb2 = 0, lbh2 = 0, lbv2 = 0;

  copyTempX(1);

  // primeira parte da heurística, pega os horizontais primeiro.
  for (int i = 1; i <= m; i++) {
    vector <int> disponiveis;
    disponiveis.clear();
    int pesoRestante = bh[i];
    int valorFixo = 0;
    for (int j = 1; j <= m; j++) {
      if (xh[i][j] + xv[i][j] != 1) {
        xh[i][j] = 0; xv[i][j] = 0;
        disponiveis.push_back(j);
      }
      else if (xh[i][j] == 1 && xv[i][j] == 0){
        valorFixo += r[i][j];
        pesoRestante -= mh[i][j];
      }
    }
    lbh += mochilaHorizontal(i, pesoRestante, disponiveis) + valorFixo;
  }

  // verificar quais os novos selecionados e fixar.
  for (int j = 1; j <= m; j++) {
    vector <int> disponiveis;
    disponiveis.clear();
    int pesoRestante = bv[j];
    int valorFixo = 0;
    for (int i = 1; i <= m; i++) {
      if (xh[i][j] == 0 && xv[i][j] == 0) { // só estão disponíveis os que ainda não foram selecionadas para a passada vertical
        disponiveis.push_back(i);
      }
      else if (xh[i][j] == 0 && xv[i][j] == 1) {
        valorFixo += r[i][j];
        pesoRestante -= mv[i][j];
      }
      else if (xh[i][j] == 1 && xv[i][j] == 1) { // se já foi selecionado pelo horizontal
        xv[i][j] = 0; // então sai desta mochila
      }
    }
    lbv += mochilaVertical(j, pesoRestante, disponiveis) + valorFixo;
  }


#ifdef DEBUG
  int totalGain = 0, tgh = 0, tgv = 0;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      if (xh[i][j] == 1)
        tgh += r[i][j];
      if (xv[i][j] == 1) 
        tgv += r[i][j];
      assert(xv[i][j] + xh[i][j] <= 1 && "Item selecionado duas vezes.");
    }
  }
  totalGain = tgh + tgv;
#endif

  assert(tgv == lbv && "Lower bound vertical não bate com o calculado.");
  assert(tgh == lbh && "Lower bound horizontal não bate com o calculado.");

  lb = lbv+lbh;
  assert(totalGain == lb);

  
  copyTempX(2); // temp_xv(h)2 tem a solucao da primeira heuristica
  copyTempXBack(1); // xv(h) tem a solucao infeasible.

  // segunda parte, pegando na vertical primeiro
  
  for (int j = 1; j <= m; j++) {
    vector <int> disponiveis;
    disponiveis.clear();
    int pesoRestante = bv[j];
    int valorFixo = 0;
    for (int i = 1; i <= m; i++) {
      if (xh[i][j] + xv[i][j] != 1) {
        xh[i][j] = 0; xv[i][j] = 0;
        disponiveis.push_back(i);
      }
      else if (xv[i][j] == 1 && xh[i][j] == 0){
        valorFixo += r[i][j];
        pesoRestante -= mv[i][j];
      }
    }
    lbv2 += mochilaVertical(j, pesoRestante, disponiveis) + valorFixo;
  }

  // verificar quais os novos selecionados e fixar.
  for (int i = 1; i <= m; i++) {
    vector <int> disponiveis;
    disponiveis.clear();
    int pesoRestante = bh[i];
    int valorFixo = 0;
    for (int j = 1; j <= m; j++) {
      if (xh[i][j] == 0 && xv[i][j] == 0) { // só estão disponíveis os que ainda não foram selecionadas para a passada horizontal
        disponiveis.push_back(j);
      }
      else if (xv[i][j] == 0 && xh[i][j] == 1) {
        valorFixo += r[i][j];
        pesoRestante -= mh[i][j];
      }
      else if (xh[i][j] == 1 && xv[i][j] == 1) { // se já foi selecionado pelo horizontal
        xh[i][j] = 0; // então sai desta mochila
      }
    }
    lbh2 += mochilaHorizontal(i, pesoRestante, disponiveis) + valorFixo;
  }

#ifdef DEBUG
  totalGain = 0, tgh = 0, tgv = 0;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      if (xh[i][j] == 1)
        tgh += r[i][j];
      if (xv[i][j] == 1) 
        tgv += r[i][j];
      assert(xv[i][j] + xh[i][j] <= 1 && "Item selecionado duas vezes.");
    }
  }
  totalGain = tgh + tgv;
#endif

  assert(tgv == lbv2 && "Lower bound vertical não bate com o calculado.");
  assert(tgh == lbh2 && "Lower bound horizontal não bate com o calculado.");

  lb2 = lbv2+lbh2;
  assert(totalGain == lb2);

  double ret;
  if (lb > lb2) {
    ret = lb;
    copyTempXBack(2); // xv e xh tem a solucao da primeira heuristica
  }
  else {
    ret = lb2;
  }
  
  return ret;
}

void updateG(int G[MAXM][MAXM], int *sumSquaredG) {
  *sumSquaredG = 0;
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      G[i][j] = 1 - xh[i][j] - xv[i][j];
      // if (G[i][j] > 0 && equal(lambda[i][j], 0.0)) {
      //   G[i][j] = 0;
      // }
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
      lambda[i][j] = min((double)(r[i][j]), max(0.0, lambda[i][j] - T*G[i][j]));
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
  double zLB_current; // uma solução viável via heurística.
  double zUP_current; // limitante dual (superior) corrente
  int G[MAXM][MAXM];
  int sumSquaredG;

  zUP = DBL_MAX; // infinito

  // aloca matriz com os lambdas
  lambda = new double*[m+1];
  for (int i = 0; i <= m; i++) {
    lambda[i] = new double[m+1]; 
  }

  initializeLambdas(); 

  zLB = 0; // ainda não temos um bom lower bound, então fica como 0
  int iteracoesSemMelhora = 0;
  // int printTimeAgain = 60;
  while (iter++ < maxiter) {

    iteracoesSemMelhora++;

    // if (!(iter % 50)) {
      // cerr << maxiter - iter << " iterações faltando." << endl;
      // cerr << "Melhor dual: " << setiosflags(ios::fixed) << setprecision(6) << zUP << endl;
      // cerr << "Melhor primal: " << (int)zLB << endl;
    // }

    // if (currentTime-startTime > printTimeAgain) {
    //   cerr << (currentTime-startTime)/60 << " minutos de execução." << endl;
    //   cerr << "Melhor dual: " << setiosflags(ios::fixed) << setprecision(6) << zUP << endl;
    //   cerr << "Melhor primal: " << (int)zLB << endl;
    //   printTimeAgain += 60;
    // }

    zUP_current = solveLR(); // resolve relaxação lagrangeana
    // cout << zUP_current << endl;
    if (zUP_current < zUP) {
      currentTime = time(NULL);
      dualTime = currentTime - startTime;
      zUP = zUP_current;
      // cerr << "Novo limitante dual: " << setiosflags(ios::fixed) << setprecision(6) << zUP << endl;
      iteracoesSemMelhora = 0;
    }
    updateG(G, &sumSquaredG); // atualiza os subgradientes e calcula a soma dos Gi^2
    T = pi*(1.05*zUP - zLB) / sumSquaredG; // atualiza o tamanho do passo
    updateLambdas(T, G);
    pi = 0.999 * pi;
    zLB_current = heuristicLB();
    if (zLB_current > zLB) {
      currentTime = time(NULL);
      primalTime = currentTime - startTime;
      zLB = zLB_current;
      copiaMatrizesX();
      // cerr << "Novo limitante primal: " << (int)zLB << endl;
    }
    assert((zLB < zUP || equal(zLB, zUP)) && "Limitante Primal maior que Dual!");

    // cout << iter << " ";
    // cout << setiosflags(ios::fixed) << setprecision(6) << zUP << " ";
    // cout << (int)zLB << endl;

    // Se tiver nSemMelhora iteracoes sem melhora, pi = pi/2;
    if (iteracoesSemMelhora == nSemMelhora) {
      // cerr << nSemMelhora << " iterações sem melhora. pi: " << pi << " -> " << pi/2 << endl; 
      pi /= 2;
      iteracoesSemMelhora = 0;
    }

    // criterios de parada
    if (equal(floor(zLB), floor(zUP))) { // solução tem que ser inteira.
      // cout << "Ótimo na iteração " << iter << endl;
      iter = maxiter;
    }
    currentTime = time(NULL);
    if (currentTime - startTime > maxTime) { // estourou o tempo
      // cout << "Tempo máximo de execução excedido: " << currentTime - startTime << " segundos" << endl;
      // cout << "Número de iterações executadas: " << iter << endl;      
      iter = maxiter;
    }
    if (pi < menorPi) {
      // cout << "Passo ficou muito pequeno." << endl;
      // cout << "Número de iterações executadas: " << iter << endl;      
      iter = maxiter;
    }    




  }

  // desaloca matriz de lambdas
  for (int i = 0; i <= m; i++) {
    delete[] lambda[i];
  }
  delete[] lambda;
}

void writeOutput(char *f) {
  ofstream output(f);
  
  output << setiosflags(ios::fixed) << setprecision(6) << zUP << endl;
  output << (int)zLB << endl;

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= m; j++) {
      if (best_xh[i][j] == 1) {
        output << i << " " << j << " h" << endl;
      }
      if (best_xv[i][j] == 1) {
        output << i << " " << j << " v" << endl;
      }
    }
  }

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
    param >> nSemMelhora; 
    param >> maxTime;
    param >> menorPi;
  }
  catch (ios_base::failure) {
    cout << "Erro na leitura do arquivo de parâmetros." << endl;
    return false;
  }

  signal(SIGINT, quit); // Liga o sinal que pára a execução sem quebrar tudo.

  startTime = time(NULL);
  currentTime = time(NULL);

  subgradientOpt(pi);
  
  writeOutput(argv[2]);

  // cout << "dualTime: " << dualTime << endl;
  // cout << "primalTime: " << primalTime << endl;

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
