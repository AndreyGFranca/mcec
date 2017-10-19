#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>
#include <string>
#include "mtwist.h"
#include <sstream>
#include <cstdio>
#include <ctime>
#include <sys/time.h>
#include <complex>

// Macro para conversão de string para numero
#define SSTR( x ) static_cast< std::ostringstream & >( \
		( std::ostringstream() << std::dec << x ) ).str()

typedef struct Disk{
	float x, y;
	unsigned int id;
} Disk;

typedef struct Celula{
	std::vector<Disk*> lista_discos;

	struct Celula* right, 
				 *left, 
				 *upper, 
				 *down, 
				 *up_right, 
				 *up_left, 
				 *down_right, 
				 *down_left;

	unsigned int linha, coluna;

}Celula;

#define PI (3.14159265358979323846)

const int 			samples = 5;
const int           Q = 30000;
const int         	N = pow(64, 2);
const int         	N_sqrt = sqrt(N) + 0.5;
const int 			vezes = 300;
const double      	eta = 0.705;
const double      	sigma = sqrt(eta / (N * PI));
const double      	D = N_sqrt;
const double      	delxy = 1.0/ (2.0 * N_sqrt);
const double        two_delxy = 2.0 * delxy;


void delx_dely(double L_jx, double L_jy, double L_kx, double L_ky, double *dx, double *dy)
{
	double d_x, d_y;

	d_x = fmod((L_jx - L_kx), 1.0);
	if (d_x > 0.5) d_x -= 1.0;
	d_y = fmod((L_jy - L_ky), 1.0);
	if (d_y > 0.5) d_y -= 1.0;

	*dx = d_x;
	*dy = d_y;
}

double dist(double L_kx, double L_ky, double L_jx, double L_jy)
{
	double d_x, d_y;

	d_x = fmod(fabs(L_kx - L_jx), 1.0);
	d_x = fmin(d_x, 1.0 - d_x);

	d_y = fmod(fabs(L_ky - L_jy), 1.0);
	d_y = fmin(d_y, 1.0 - d_y);


	return sqrt(fabs(pow(d_x, 2) + pow(d_y, 2)) );
}

struct DistanceFunc
{
	DistanceFunc(const Disk* _d) : d(_d) {}

	bool operator()(const Disk* lhs, const Disk* rhs) const
	{
		return dist(d->x, d->y, lhs->x, lhs->y) < dist(d->x, d->y, rhs->x, rhs->y);
	}

	private:
	const Disk* d;
};


std::complex<double> calc_psi_k(Disk* disk_k, 
	Celula cell_k,
	unsigned int k)
{
	int n_neighbor = 0;
	float gamma = 2.8;
	std::vector<Disk*> k_list;
	double dx, dy, angle;
	Celula* aux;
	std::complex<double> vetor(0.0, 0.0);

	for (int j = 0; j < 17; ++j)
	{
		if(j == 0)
			aux = cell_k.right;
		else if(j == 1)
			aux = cell_k.up_right;
		else if (j == 2)
			aux = cell_k.upper;
		else if (j == 3)
			aux = cell_k.up_left;
		else if(j == 4)
			aux = cell_k.left;
		else if(j == 5)
			aux = cell_k.down_left;
		else if (j == 6)
			aux = cell_k.down;
		else if (j == 7)
			aux = cell_k.down_right;
		else if (j == 8)
			aux = &cell_k;
		if(j == 9)
			aux = cell_k.right->right;
		else if(j == 10)
			aux = cell_k.up_right->up_right;
		else if (j == 11)
			aux = cell_k.upper->upper;
		else if (j == 12)
			aux = cell_k.up_left->up_left;
		else if(j == 13)
			aux = cell_k.left->left;
		else if(j == 14)
			aux = cell_k.down_left->down_left;
		else if (j == 15)
			aux = cell_k.down->down;
		else if (j == 16)
			aux = cell_k.down_right->down_right;

		for(int i = 0; i < aux->lista_discos.size(); i++){
			k_list.push_back(aux->lista_discos[i]);
		}
	}

	std::sort(k_list.begin(), k_list.end(), DistanceFunc(disk_k));

	for (int i = 0; i < k; ++i)
	{
		if(disk_k->id !=  k_list[i]->id){
			n_neighbor++;
			delx_dely(k_list[i]->x, k_list[i]->y, disk_k->x, disk_k->y, &dx, &dy);
			angle = std::arg(std::complex<double>(dx, dy));
			vetor += std::exp(std::complex<double>(0.0,6.0) * std::complex<double>(angle, 0));
		}
	}

	if (n_neighbor > 0)
	{
		vetor /= (double)n_neighbor;
		//std::cout << "Numero de Vizinhos = " << n_neighbor << std::endl;
		///std::cout << "Vetor  = " << vetor << std::endl;
	}
	return vetor;
}

std::complex<double> calc_psi_global(Disk (&disk)[N_sqrt][N_sqrt], Celula (&celula)[N_sqrt][N_sqrt])
{
	std::complex<double> sum_vetor(0.0, 0.0);
	for (int i = 0; i < N_sqrt; ++i)
	{
		for (int j = 0; j < N_sqrt; ++j)
		{
			unsigned int ci = ceil(disk[i][j].x / two_delxy) - 1; // coluna
			unsigned int cj = ceil(disk[i][j].y / two_delxy) - 1;
			sum_vetor += calc_psi_k(&disk[i][j], celula[ci][cj], 6);
		}
	}

	return sum_vetor / (double)N;
}

void novoL(Disk (&disk)[N_sqrt][N_sqrt], 
	Celula (&celula)[N_sqrt][N_sqrt],
	double D, 
	double sigma,
	int N_sqrt);

double event(double b_x, 
	double b_y, 
	double a_x, 
	double a_y, 
	int dirc, 
	double sigma);


time_t timer;
char buffer[26];
struct tm* tm_info;

int main(){

	double Psi_mod[Q+1], Psi_mod_sq[Q+1], erroPsi[Q+1];
	std::complex<double> Psi(0.0, 0.0);

	char nome4[60];
	sprintf(nome4, "resultados/resultados_m%d_eta%.f_Q%d_S%d.csv", N_sqrt,1000*eta,Q,samples);
	FILE *resultados   = fopen(nome4,"w");

	struct timeval start, end;
	gettimeofday(&start, NULL);

	Celula 			celula[N_sqrt][N_sqrt];
	Disk 			disk[N_sqrt][N_sqrt];
	//mt_state 		state[50];
	double duracao;

	time(&timer);
	tm_info = localtime(&timer);

	strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
	puts(buffer);

	unsigned int id_count = 1;

	printf("Raiz do numero de discos N0: %d\n", N_sqrt);
	printf("Numero de iteracoes Q: %d\n", Q);
	printf("Densidade eta: %f\n", eta);
	printf("Numero de amostras: %d\n", samples);
	printf("Numero de vezes: %d\n", vezes);

	/****************************************************
	 * 			  DESCOBRINDO AS VIZINHAS 				*
	 * Atribuindo a cada célula o endereço de suas célul*
	 * as vizinhas.										*
	 ****************************************************/
	for (int i = 0; i < N_sqrt; ++i)
	{
		for (int j = 0; j < N_sqrt; ++j)
		{
			// Definição dos vizinhos

			// Se for a primeira  célula da grade
			if(i == 0 && j == 0){
				// Armazenando a linha e a coluna da celula
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][N_sqrt-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[N_sqrt-1][j];
				celula[i][j].down_right = &celula[N_sqrt-1][j+1];
				celula[i][j].down_left = &celula[N_sqrt-1][N_sqrt-1];
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][N_sqrt-1];
			}
			else if(i == 0 && j == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[0][0];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[N_sqrt-1][j];
				celula[i][j].down_right = &celula[N_sqrt-1][0];
				celula[i][j].down_left = &celula[N_sqrt-1][j-1];
				celula[i][j].up_right = &celula[i+1][0];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
			else if(i==N_sqrt-1 && j == 0){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][N_sqrt-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[0][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][N_sqrt-1];
				celula[i][j].up_right = &celula[0][j+1];
				celula[i][j].up_left = &celula[N_sqrt-1][N_sqrt-1];
			}
			else if(i == N_sqrt-1 && j == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][0];
				celula[i][j].upper = &celula[0][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][0];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[0][0];
				celula[i][j].up_left = &celula[0][j-1];
			}
			else if(i > 0 && i <= N_sqrt-2 && j == 0){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][N_sqrt-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][N_sqrt-1];
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][N_sqrt-1];
			}
			else if(j > 0 && j < N_sqrt-1 && i == 0){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[N_sqrt-1][j];
				celula[i][j].down_right = &celula[N_sqrt-1][j+1];
				celula[i][j].down_left = &celula[N_sqrt-1][j-1]; 
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
			else if(i > 0 && i <= N_sqrt-2 && j == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][0];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][0];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[i+1][0];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
			else if(j > 0 && j <= N_sqrt-2 && i == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[0][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[0][j+1];
				celula[i][j].up_left = &celula[0][j-1];
			}
			else{
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
		}
	}


	/****************************************************
	 * 			  ADCIONANDO DISCOS AS CELULAS  			*
	 * Adciona a cada celula, o seu respectivo disco 	*
	 ****************************************************/
	for (int i = 0; i < N_sqrt; ++i){
		for (int j = 0; j < N_sqrt; ++j){
			celula[i][j].lista_discos.push_back(&disk[j][i]);
		}
	}

	//gerando configuração inicial quadrado
	// for (int i = 0; i < N_sqrt; i++){
	// 	for (int j = 0; j < N_sqrt; j++){
	// 		disk[i][j].x = delxy + i * two_delxy;
	// 		disk[i][j].y = delxy + j * two_delxy;
	// 		disk[i][j].id = id_count++;
	// 	}
	// }

	// Gravando a configuracao inicial em arquivo
	//std::ofstream myfile;
	//std::string s1 = "resultados/LxLy_inicial_N0" + SSTR( N_sqrt ) + "_eta" + SSTR( 100*eta ) + "_se"+SSTR(ks) + "_Q" + SSTR( Q ) + ".csv";
	//myfile.open (s1.c_str());
	//for (int i = 0; i < N_sqrt; ++i)
	//{
		//for (int j = 0; j < N_sqrt; ++j)
		//{
			//myfile << disk[i][j].x << "," << disk[i][j].y << std::endl;
		//}
	//}
	//myfile.close();

	//gerando configuração inicial triangular
	// int cont = 0;
	// for (int i = 0; i < N_sqrt; i++){
	// 	for (int j = 0; j < N_sqrt; j++){
	// 		if(j%2 == 0){
	// 			disk[i][j].x = delxy + i * two_delxy;
	// 			disk[i][j].y = delxy + j * two_delxy;
	// 			disk[i][j].id = id_count++;
	// 		}
	// 		else{
	// 			disk[i][j].x = (delxy + i * two_delxy) + delxy;
	// 			disk[i][j].y = delxy + j * two_delxy;
	// 			disk[i][j].id = id_count++;
	// 		}
	// 	}
	// }

	for (int k = 0; k < samples; k++){

		mt_seed32(k);
		// Limpa a lista de discos
		for (int i = 0; i < N_sqrt; ++i){
			for (int j = 0; j < N_sqrt; ++j){
				celula[i][j].lista_discos.clear();
			}
		}

		//adciona
		for (int i = 0; i < N_sqrt; ++i){
			for (int j = 0; j < N_sqrt; ++j){
					celula[i][j].lista_discos.push_back(&disk[j][i]);
			}
		}

		// gera configuracao inicial quadrada
		// for (int i = 0; i < N_sqrt; i++){
		//     for (int j = 0; j < N_sqrt; j++){
		//         disk[i][j].x = delxy + i * two_delxy;
		//         disk[i][j].y = delxy + j * two_delxy;
		//         disk[i][j].id = id_count++;
		//     }
		// }


		//gerando configuração inicial triangular
		int cont = 0;
		for (int i = 0; i < N_sqrt; i++){
			for (int j = 0; j < N_sqrt; j++){
				if(j%2 == 0){
					disk[i][j].x = delxy + i * two_delxy;
					disk[i][j].y = delxy + j * two_delxy;
					disk[i][j].id = id_count++;
				}
				else{
					disk[i][j].x = (delxy + i * two_delxy) + delxy;
					disk[i][j].y = delxy + j * two_delxy;
					disk[i][j].id = id_count++;
				}
			}
		}

		//std::cout << "iteracao = " << k << std::endl;
		//printf("Iniciando amostra: %d\n", k);
		for (int i = 0; i <= Q; i++){
			if(i % vezes == 0){
				printf("Iniciando iteracao: %d\n", i);
				//std::ofstream myfile;
				//std::string s = SSTR( i );
				//std::string s2 = "resultados/LxLy_N0" + SSTR( N_sqrt )  + "_eta" + SSTR( 100*eta )+ "_t" + SSTR( i )  + ".csv";
				//myfile.open (s2.c_str());
				//for (int i = 0; i < N_sqrt; ++i){
					//for (int j = 0; j < N_sqrt; ++j){
						//myfile << disk[i][j].x << "," << disk[i][j].y << std::endl;
					//}
				//}
				//myfile.close();
				Psi = calc_psi_global(disk, celula);
				Psi_mod[i] +=  std::abs(Psi) / (double) samples;
				Psi_mod_sq[i] +=  std::pow( Psi_mod[i], 2) / (double) samples; // SAMPLES;
			}
			novoL(disk, celula, D, sigma, N_sqrt);
		} // fecha o for do tempo
	}

// calculo do erro em Psi
	for(int i = 0; i <= Q; i++){
		erroPsi[i] = std::sqrt(std::fabs(Psi_mod[i] - std::pow(Psi_mod_sq[i], 2))); /// (SAMPLES - 1)));
	}
	
//gravando os valores de Psi em um arquivo
	for (int i4 = 0; i4 <= Q; i4++){
		if (i4 % vezes == 0){
			fprintf(resultados, "%d, %f, %f\n", i4, Psi_mod[i4], erroPsi[i4]);
		}
	}
	////Gravando configuracao final em um arquivo
	//std::ofstream myfile2;
	//std::string s3 = "resultados/LxLy_final_N0" + SSTR( N_sqrt ) + "_eta" + SSTR( 100*eta ) + "_se"+SSTR(ks) + "_Q" + SSTR( Q ) + ".csv";
	//myfile2.open (s3.c_str());
	//for (int i = 0; i < N_sqrt; ++i)
	//{
		//for (int j = 0; j < N_sqrt; ++j)
		//{
			//myfile2 << disk[i][j].x << "," << disk[i][j].y << std::endl;
		//}
	//}
	//myfile2.close();

	//Tempo final do programa.
	gettimeofday(&end, NULL);
	duracao = (((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6);

	printf("Tempo de execucao em segundos: %f\n", duracao);
	printf("Tempo de execucao em horas: %f\n", duracao / 3600);
	printf("Tempo de execucao por iteracao e por amostra em segundos: %f\n", duracao / (samples*Q));

	char nome3[60];

	sprintf(nome3, "resultados/parametros_m%d_eta%.f_Q%d_S%d.txt", N_sqrt,1000*eta,Q,samples);
	FILE *parametros   = fopen(nome3,"w");

	fprintf(parametros, "Parametros importantes.\n");
	fprintf(parametros, "Numero de discos N: %d\n", N);
	fprintf(parametros, "Raiz de N = m: %d\n", N_sqrt);
	fprintf(parametros, "Numero de iteracoes Q: %d\n", Q);
	fprintf(parametros, "Densidade eta: %f\n", eta);
	fprintf(parametros, "Numero de vezes: %d\n", vezes);
	fprintf(parametros, "Numero de amostras: %d\n", samples);
	fprintf(parametros, "Tempo de execucao: %f\n", duracao);
	fprintf(parametros, "Tempo de execucao em horas: %f\n", duracao / 3600);
	fprintf(parametros, "Tempo de execucao por iteracao em segundos: %f\n", duracao / (samples*Q));

	//sprintf(nome4, "resultados/valores.csv");
	//FILE *valores   = fopen(nome4,"w");
	//fprintf(valores, "%d,%f,%d,%d,%d\n", N_sqrt, eta, Q,vezes, ks);
	
	printf("Fim.\n");

}

void novoL(Disk (&disk)[N_sqrt][N_sqrt], Celula (&celula)[N_sqrt][N_sqrt], double D, double  sigma, int N_sqrt)
{
	double event_min = 0;

	 //Gerando direcao aleatoria 0 ou 1

	int dirc = (unsigned int) mt_lrand() % 2;

	double distance_to_go = D;
	Disk *next_a;

	int rand_i1, rand_i2;

	 //Gerando posicoes aleatorias.
	rand_i1 = (unsigned int) mt_lrand() % N_sqrt ;
	rand_i2 = (unsigned int) mt_lrand() % N_sqrt ;

	next_a = &disk[rand_i1][rand_i2];
	while(distance_to_go > 0.0)
	{
		bool flag = false;
		Disk* a = next_a;
		event_min = distance_to_go;

		// Descobrimos em qual linha e coluna o disco está:
		unsigned int l = ceil(a->y / two_delxy) - 1; // linha
		unsigned int c = ceil(a->x / two_delxy) - 1; // coluna

		
		// Calcula a distancia para todos os discos das celulas vizinhas na direcao em que 
		{
		Celula* aux;

		for (int i = 0; i < celula[l][c].lista_discos.size(); i++){
			// Se o disco a verificar não for o disco a:
			if(celula[l][c].lista_discos[i]->id != a->id){
				double event_b = event(celula[l][c].lista_discos[i]->x, 
					celula[l][c].lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = celula[l][c].lista_discos[i];
				}

			}
		}
		// Celula da Direita
		aux =  celula[l][c].right;
		for (int i = 0; i < 3; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->right;
		}

		// Celula upper-right
		aux =  celula[l][c].up_right;
		for (int i = 0; i < 3; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = aux->lista_discos[i];
				}
			}
			if(dirc == 0){
				aux = aux->upper;
			}
			else{
				aux = aux->right;
			}
		}

		aux =  celula[l][c].down_right;
		for (int i = 0; i < 2; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->right;
		}

		aux =  celula[l][c].upper;
		for (int i = 0; i < 2; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->upper;
		}
		aux = celula[l][c].down;
		for (int i = 0; i < 2; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->down;
		}

		aux = celula[l][c].up_left;
		for (int i = 0; i < 2; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
//					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->upper;
		}

		aux = celula[l][c].left;
		for (int i = 0; i < 2; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
//					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->left;
		}

		aux = celula[l][c].down_left;
		for (int i = 0; i < 2; ++i)
		{
			for (int i = 0; i < aux->lista_discos.size(); i++){
				double event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					next_a = aux->lista_discos[i];
				}
			}
			aux = aux->down;
		}


		}


		if (dirc == 1){
			//Verifica qual celula o disco estava;
			unsigned int l_anterior = l;
			unsigned int c_anterior = c;
			//Atualiza a posicao
			a->x = std::fmod((a->x + fabs(event_min-0.0001)), 1.0);
			unsigned int c_atual = ceil(a->x / two_delxy) - 1; // coluna
			unsigned int l_atual = ceil(a->y / two_delxy) - 1;
			if(c_anterior != c_atual || l_atual != l_anterior){
				Disk* d;
				//std::cout << "O Disco atualizou sua posição e mudou de célula:\nCelula Anterior: ("<<l<<", "<<c<<")"<<std::endl;
				//std::cout << "Celula Atual: (" <<l_atual << ", "<<c_atual<<")"<<std::endl; 
				for (int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
					if(celula[l_anterior][c_anterior].lista_discos[i]->id == a->id){
						d = celula[l_anterior][c_anterior].lista_discos[i];
						celula[l_anterior][c_anterior].lista_discos.erase(celula[l_anterior][c_anterior].lista_discos.begin() + i);
						break;
					}
				}
				//std::cout << "Verificando se removeu da celula anterior, e atualizou:" << std::endl;
				for(int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
				}

				//Vou na célula atual e armazeno o disco:
				celula[l_atual][c_atual].lista_discos.push_back(d);
				for(int i = 0; i < celula[l_atual][c_atual].lista_discos.size(); ++i){
				}
				
			}
		}
		else{
			unsigned int l_anterior = l;
			unsigned int c_anterior = c;
			a->y = std::fmod((a->y + fabs(event_min-0.0001)), 1.0);
			//Verifica novamente a celula
			unsigned int c_atual = ceil(a->x / two_delxy) - 1; // coluna
			unsigned int l_atual = ceil(a->y / two_delxy) - 1;
			//Se ele mudou de celula
			if(l_anterior != l_atual || l_atual != l_anterior){
				Disk* d;
				for (int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
					if(celula[l_anterior][c_anterior].lista_discos[i]->id == a->id){
						d = celula[l_anterior][c_anterior].lista_discos[i];
						celula[l_anterior][c_anterior].lista_discos.erase(celula[l_anterior][c_anterior].lista_discos.begin() + i);
						break;
					}
				}
				//std::cout << "Verificando se removeu da celula anterior, e atualizou:" << std::endl;
				//std::cout << "	 ==========================================	 " << std::endl;
				for(int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
				}

				//Vou na célula atual e armazeno o disco:
				celula[l_atual][c_atual].lista_discos.push_back(d);
				for(int i = 0; i < celula[l_atual][c_atual].lista_discos.size(); ++i){
				}
				//std::cout << "	 ==========================================	 " << std::endl;
			}
		}
		distance_to_go -= event_min;

	}


}


double event(double b_x, double b_y, double a_x, double a_y, int dirc, double sigma)
{
	double d_perp = 0.0;
	double d_para = 0.0;

	if (dirc == 1)
		d_perp = std::fmod(fabs(b_y - a_y), 1.0);
	else
		d_perp = std::fmod(fabs(b_x - a_x), 1.0);

	d_perp = fmin(d_perp, 1.0 - d_perp);
	if (d_perp > 2.0 * sigma)
		return (double)INFINITY;
	else{
		d_para = sqrt(fabs(4.0 * sigma*sigma - d_perp*d_perp));

		if (dirc == 1){
			return std::fmod((b_x - a_x - d_para + 1.0), 1.0);
		}
		else{
			return std::fmod((b_y - a_y - d_para + 1.0), 1.0);
		}
	}
}
