#include <iostream>
#include <fstream>
#include <vector>

#define MAXM 1000
#define MAXV 1000

using namespace std;

int m, k;
int bh[MAXM], bv[MAXM]; // capacidade de memória do satélite para as faixas
int r[MAXM][MAXM]; // ganho pela captura do fragmento i j
int mh[MAXM][MAXM], mv[MAXM][MAXM]; // memória necessária para capturar fragmento ij em cada observacao.
int v; // numero de vertices do poligono
int vx[MAXV], vy[MAXV];



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
    cout << p+1<< " " << Fh[p+1] << " ";
  }
  cout << endl;
  for (int p = 0; p < m; p++) { 
    cout << p+1 << " " << Fv[p+1] << " ";
  }
  cout << endl << k << endl;
}

// // Resolve a mochila da faixa i
// void solveLRi(int i) {


// }

int main(int argc, char **argv) {

  if (argc < 3) {
    cout << "USAGE: " << argv[0] << " ENTRADA SAIDA" << endl;
  }

  // tenta ler entrada
  bool r;
  r = readInput(argv[1]);
  if (!r) {
    cout << "Erro ao ler entrada!" << endl;
    return 1;
  }
  
  

  return 0;
}
