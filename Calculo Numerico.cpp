#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include <iostream>

using namespace std;


#define N 2 //Tamanho da matriz pro gauss e parcial



//VC define as expressoes antes de execultar o codigo pra responder o estilo gauss-seider e jacobi
void definesUsados(){
	
	//Define o sistema linear que ir� ser resolvido Gauss-Seider e Jacobi
	#define f1(x,y,z)  (17-y+2*z)/20
	#define f2(x,y,z)  (-18-3*x+z)/20
	#define f3(x,y,z)  (25-2*x+3*y)/20
	
	//Define os valores da linha 0 da primeira itera��o
	#define aux1 1
	#define aux2 1.5
	#define aux3 2
}





//Imprime Matriz de Gauss
void imprimeMatriz(float M[N][N+1]){
	cout << "\n";
	for(int i = 0; i<N; i++){
        for(int j = 0; j<=N; j++){
            printf("|%.3f\t", M[i][j]);
        }
      	cout << "\n";
    }
}

//Imprime matriz para pivotamento parcial
void imprimeMatriz2(float M[N][N], float X[N]){
	cout << "\n";
	
	//la�o que percorre a matriz e o vetor e imprime na tela
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("|%.3f ", M[i][j]);
		}
		printf(" |%.3f", X[i]);
		cout << "\n";
	}
}






//Fun��o que percorre a matriz e troca as linhas da mesma 
void trocaLinhas(int l1, int l2, float M[N][N], float X[N]){
		double aux = 0;
		
		// Troca elementos das linhas na matriz M[][]
		for (int i = 0; i < N; i++) {
			aux = M[l1][i];
			M[l1][i] = M[l2][i];
			M[l2][i] = aux;
		}
		
		//Troca os elementos contidos no vetor X[][]
		aux = X[l1];
		X[l1] = X[l2];
		X[l2] = aux;

	}
	
	
	
	
	
	
//Funcao para calcular matriz usando metodo de Pivotamento Parcial	
void pivotamentoParcial(float M[N][N], float X[N]){
		int k, i;
		// Percorre as linhas o la�o se repete at� n�mero de linhas da matriz - 1
		for (k = 0; k < (N - 1); k++) {
			// Percorre linhas secund�rias at� N
			for (i = (k + 1); i < N; i++) {
				// Escolhe o piv� do elemento de maior m�dulo entre os elementos da matriz
				//OBS: fun��o fabs() faz o m�dulo de um n�mero float
				if (fabs(M[k][k]) < fabs(M[i][k])) {
					//chama a fun��o de troca de linhas
					trocaLinhas(k, i, M, X); 

				}
				
				//elimina��o de gauss

				double fator = M[i][k] / M[k][k];
				M[i][k] = 0;

				// Calcula as pr�xima linhas com base nos fatores de multiplica��o encontrados...
				int j;
				for (j = (k + 1); j < N; j++) {
					M[i][j] = M[i][j] - fator * M[k][j];
				}

				// C�lculo do vetor resutado
				X[i] = X[i] - fator * X[k];
			}
		}
	}

//Funcao para calcular matriz usando metodo de Gauss	
void eliminacaoProgressiva(float M[N][N+1]){
	float fator; 
	
    for(int k = 0; k<N-1; k++){
	//Coluna pivo. 
        for(int i = (k+1); i<N; i++){
		//Linhas que sofrerao as eliminacoes 
		//Sempre iniciarao na linha seguinte da linha pivo
            fator = M[i][k] / M[k][k];
			//fator multiplicador
         
		    for(int j=0; j<=N; j++){
			//Coluna da linha que sofrerao modificacoes da linha i
		        M[i][j] = M[i][j] - fator*M[k][j];
            }
        }
    }
}

//Funcao para calcular matriz usando metodo de Gauss-Seider
void gaussSeider(){
	
	definesUsados();
	
	float x0 = aux1, y0 = aux2, z0 = aux3, x1, y1, z1, e1, e2, e3, e;  //valor das variaveis na primeira itera��o
	 int cont = 1; //inicializando a vari�vel contadora 
	 
	 cout << "Digite o valor do erro: " << endl; //Leitura do valor do erro necess�rio no crit�rio de parada 
	 cin >> e;
	 
	 cout << "\nIteracao\tx\ty\tz" << endl;
	 
	 //inicio do loop de itera��es
	 printf("\n%d\t%0.4f\t%0.4f\t%0.4f\n",0, x0,y0,z0); 
	 
	 do
	 {
	  // realiza o c�lculo da primeira itera��o atribuindo 0
	  x1 = f1(x0,y0,z0);
	  y1 = f2(x1,y0,z0);
	  z1 = f3(x1,y1,z0);
	  
	  printf("%d\t%0.4f\t%0.4f\t%0.4f\n",cont, x1,y1,z1);
	
	  //realiza o c�lculo do erro
	  e1 = fabs(x0-x1);
	  e2 = fabs(y0-y1);
	  e3 = fabs(z0-z1);
	
	  cont++;		//vari�vel contadora de itera��es
	
	  //seta o valor das pr�ximas itera��es 
	  x0 = x1;
	  y0 = y1;
	  z0 = z1;
	
	 }while(e1>e && e2>e && e3>e); // teste para saber se atingiu o crit�rio de parada
	 
	 //imprime a solu��o
	 printf("\nSolucao: x=%0.3f, y=%0.3f e z = %0.3f\n",x1,y1,z1);
	
}

//Funcao para calcular matriz usando metodo de Jacobi
void jacobi(){
	
	definesUsados();
	
	float x0=aux1, y0=aux2, z0=aux3, x1, y1, z1, e1, e2, e3, e; //inicializando as vari�veis 
	 int cont=1; //inicializando a vari�vel contadora das itera��es 
	 
	 //leitura do valor do erro para crit�rio de parada
	 cout << "Digite o valor do erro: ";
	 cin >> e;
	
	 
	 cout << "\nIteracao\tx\ty\tz" << endl; //cabe�alho
	 printf("\n%d\t%0.4f\t%0.4f\t%0.4f\n",0, x0,y0,z0); // print da linha o da eq. que pode ser definida no define. 
	 
	 //in�cio do loop de itera��es 
	 do
	 {
	  //realiza o c�lculo da primeira itera��o atribuindo 0 (a depender da quest�o)
	  x1 = f1(x0,y0,z0);
	  y1 = f2(x0,y0,z0);
	  z1 = f3(x0,y0,z0);
	  printf("%d\t%0.4f\t%0.4f\t%0.4f\n",cont, x1,y1,z1);
	
	  //realiza o c�lculo do erro
	  e1 = fabs(x0-x1);
	  e2 = fabs(y0-y1);
	  e3 = fabs(z0-z1);
	
	  cont++; //incremento da vari�vel contadora 
	
	  //define os valores da pr�xima itera��o
	  x0 = x1;
	  y0 = y1;
	  z0 = z1;
	  
	 }while(e1>e && e2>e && e3>e); //teste para saber se atingiu o crit�rio de parada 
	
	 printf("\nSolucao: x=%0.3f, y=%0.3f z = %0.3f\n",x1,y1,z1); //imprime a solu��o na tela
}









int main (){
    
    
    //Para Metodo de gauss passando os resultados na matriz
    /*float M[N][N+1] = {{2, 2, 1, 1, 7},
                       {1, -1, 2, -1, 1},
                       {3, 2, -3, -2, 4},
                       {4, 3, 2, 1, 12}};
					   */
                 
                 
	//Para Pivotamento parcial separando os resultados em um vetor    
	/*float M[N][N] ={{0.004,15.73},
 	 	 	 	 {0.423,-24.72}};*/
 	 	 	 	 
	//vetor dos resultados para pivotamemto parcial			
	//float X[N] = {15.77,-20.49};




	//cout << "\n\t\tMatriz original" << endl;
	//imprimeMatriz(M);
	//imprimeMatriz2(M, X);
	
	//eliminacaoProgressiva(M);
	//pivotamentoParcial(M, X);
	//gaussSeider();
	jacobi();
	
	//cout << "\n\t\tMatriz de Gauss" << endl;
    //imprimeMatriz(M);
    //cout << "\n\t\tMatriz de Pivotamento" << endl;
	//imprimeMatriz2(M, X);
    
    return 0;
}

