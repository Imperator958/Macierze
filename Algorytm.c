#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define N 4000
#define P 4
#define PP 2

float A[N][N], B[N][N], C_glob[N][N], C_sek[N][N];
float a[N / PP][N / PP], b[N / PP][N / PP], c[N / PP][N / PP];
float aa[N / PP][N / PP], bb[N / PP][N / PP];
float (*psa)[N / PP], (*psb)[N / PP], (*pra)[N / PP], (*prb)[N / PP];
float tmpa[N / PP][N / PP], tmpb[N / PP][N / PP], tmpc[N / PP][N / PP];

double startwtime1, startwtime2, endwtime;

int modulo(int a, int b){
	int r = a % b;
	return r < 0 ? r + b : r;
}

int to_right(int rank, int p, int pp, int times)
{
    if (times == 0){
        return rank;
    }
    
    int score = rank;
    
    for (int i = 0; i < times; i++){
        score = ((score + 1) % pp) + ((score / pp) * pp);
    }
    
    return score;
}

int to_down(int rank, int p, int pp, int times)
{
    if (times == 0){
        return rank;
    }
    
    int score = rank;
    
    for (int i = 0; i < times; i++){
        score = (score + PP) % P;
    }
    
    return score;
}


int main(int argc, char** argv)
{

	FILE* plik;
	FILE* plik_out;

	int my_rank, size;
	int row, col, mod = 0;
	int data_received = -1;
	int tag = 101;
	int koniec;
	int x,y,z,t;
	

	MPI_Status statSend[2], statRecv[2], statRecvCollect[P];
    	MPI_Request reqSend[2], reqRecv[2], reqSendCollect[P], reqRecvCollect[P];
    	
    
	MPI_Init(0, 0);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	row = my_rank / PP; 
	col = my_rank % PP;
	
	int gór = (my_rank + P - PP) % P;
	int dol = (my_rank + PP) % P;
	int praw = ((my_rank + 1) % PP) + ((my_rank / PP) * PP);
	int lew = modulo(my_rank - 1, PP) + ((my_rank / PP) * PP);

	if (my_rank == 0)
		printf("obliczenia metod  Cannona dla tablicy %d x %d element w \n", N, N);

	if (my_rank == 0) startwtime1 = MPI_Wtime();//czas w sekundach

	//wczytanie danych przez proces rank=0
	if (my_rank == 0)
	{
		plik = fopen("liczby.txt", "r");
		if (plik == NULL)
		{
			printf("Blad otwarcia pliku \"liczby.txt\"\n");
			koniec = 1;
			MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize();
			exit(0);
		}
		else {
			koniec = 0;
			MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (koniec) { 
			MPI_Finalize(); 
			exit(0); 
		}
	}


	if (size != P) {
		if (my_rank == 0) 
			printf("wywolano obliczenia iloczynu macierzy metoda cannona na %d procesach - uruchom mpiexec -n %d matrixmult\n", size, P);
			MPI_Finalize(); 	
			exit(0);
	}

	if (my_rank == 0)
	{
		for (int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++){
        			fscanf(plik, "%f", &A[i][j]);
        		}
  		}	
  		
  		//rewind(plik);
  		
  		for (int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++){
        			fscanf(plik, "%f", &B[i][j]);
        		}
  		}	
  		
  		for(int i = 0; i < N/PP; i++){
  			for(int j = 0 ; j < N/PP; j++)
  			{
  				a[i][j] = A[i][j];
  				b[i][j] = B[i][j];
  			}
  		
  		}
  		
  		for(int i = 1; i < P; i++)
  		{
  			for(int j = 0; j < N/PP; j++)
  			{
  				for(int k = 0; k < N/PP; k++)
  				{
					/*x = (to_right(i, P, PP, i / PP) / PP) * (n / PP) + j;
					z = (to_right(i, P, PP, i / PP) % PP) * (n / PP) + k;
					y = (to_down(i, P, PP, i % PP) / PP) * (n / PP) + j;
					t = (to_down(i, P, PP, i % PP) % PP) * (n / PP) + k;
					
					tmpa[j][k] = A[(to_right(i, P, PP, i / PP) / PP) * (n / PP) + j][(to_right(i, P, PP, i / PP) % PP) * (n / PP) + k];
					
					tmpb[j][k] = B[(to_down(i, P, PP, i % PP) / PP) * (n / PP) + j][(to_down(i, P, PP, i % PP) % PP) * (n / PP) + k];*/
					
					tmpa[j][k] = A[(to_right(i, P, PP, i / PP) / PP) * (N / PP) + j][(to_right(i, P, PP, i / PP) % PP) * (N / PP) + k];
					
					tmpb[j][k] = B[(to_down(i, P, PP, i % PP) / PP) * (N / PP) + j][(to_down(i, P, PP, i % PP) % PP) * (N / PP) + k];
		
				}
			}
			
			MPI_Isend(tmpa, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqSend[0]);
			//test konca komunikacji
			MPI_Isend(tmpb, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqSend[1]);
			//test konca komunikacji
		}

	}
	else
	{
		MPI_Irecv(a, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[0]);
		MPI_Irecv(b, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[1]);
		
		MPI_Waitall(2, reqRecv, statRecv);
		
		//test konca komunikacji
	}
	
	pra = aa; 
	prb = bb; 
	psa = a; 
	psb = b;
	//przygotowanie lokalnej  tablicy wynikowej 
	
	for (int i = 0; i < N / PP; i++)
		for (int j = 0; j < N / PP; j++)
		{
			c[i][j] = 0;
		}


	if (my_rank == 0) 
		startwtime2 = MPI_Wtime();//czas w sekundach
	
	for (int kk = 0; kk < PP; kk++)
	{
		//krok obliczeń
		for (int i = 0; i < N / PP; i++)
			for (int j = 0; j < N / PP; j++)
				for (int k = 0; k < N / PP; k++)
					c[i][j] += psa[i][k] * psb[k][j];
					
		//komunikacja pra = aa; prb = bb; psa = a; psb = b;
		MPI_Irecv(pra, N*N / PP / PP, MPI_FLOAT, praw, tag, MPI_COMM_WORLD, &reqRecv[0]);
		MPI_Irecv(prb, N*N / PP / PP, MPI_FLOAT, dol, tag, MPI_COMM_WORLD, &reqRecv[1]);		
		MPI_Isend(psa, N*N / PP / PP, MPI_FLOAT, lew, tag, MPI_COMM_WORLD, &reqSend[0]);
		MPI_Isend(psb, N*N / PP / PP, MPI_FLOAT, gór, tag, MPI_COMM_WORLD, &reqSend[1]);
		
		MPI_Waitall(2, reqSend, statSend);
        	MPI_Waitall(2, reqRecv, statRecv);
		
		if (mod = ((mod + 1) % 2))
		{ 
			pra = a; 
			prb = b; 
			psa = aa; 
			psb = bb; 
		}
		else
		{ 
			pra = aa; 
			prb = bb; 
			psa = a; 
			psb = b; 
		}
		
	//obliczenia iloczynu macierzy zgodnie z algorytmem Cannona 
	//do uzupenienia
	}
	
	if (my_rank == 0)
	{
		
		endwtime = MPI_Wtime();
		printf("Calkowity czas przetwarzania wynosi %f sekund\n", endwtime - startwtime1);
		printf("Calkowity czas obliczen wynosi %f sekund\n", endwtime - startwtime2);
	

	}
	// test poprawnosci wyniku 
	
	MPI_Isend(c, N * N / PP / PP, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqSendCollect[my_rank]);
	

	if (my_rank == 0)
	{
		// odbiór wyników obliczeń równoległych do globalnej tablicy wynikowej Cglob
		for (int i = 0; i < P; i++)
        	{
          		MPI_Irecv(tmpc, N * N / PP / PP, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqRecvCollect[i]);
      
			MPI_Wait(&reqRecvCollect[i], &statRecvCollect[i]);
			
			for (int j = 0; j < N / PP; j++){
				for (int k = 0; k < N / PP; k++){
					//x = (i / PP) * (n / PP) + j;
					//y = (i % PP) * (n / PP) + k;
                   			C_glob[(i / PP) * (N / PP) + j][(i % PP) * (N / PP) + k] = tmpc[j][k];
                		}
			}
		}
		
		// obliczenia sekwencyjne mnożenia tablic CSek=A*B
		startwtime2 = MPI_Wtime();

		for(int i = 0 ; i < N; i++){
			for(int j = 0; j < N; j++){
				C_sek[i][j] = 0;
				for(int k = 0; k < N; k++){
					C_sek[i][j] += A[i][k] * B[k][j]; 
				}
			}
		}
		
		endwtime = MPI_Wtime();
		printf("Calkowity czas obliczen sekwencyjnych wynosi %f sekund\n", endwtime - startwtime1);
			
		// porównanie poprawności obliczeń (Csek, Cglob) przy uwzględniniu progu poprawności 
		int Err = 0;
       		for (int i = 0; i < N; i++)
        	{
            		for (int j = 0; j < N; j++)
            		{
                		if (C_sek[i][j] != C_glob[i][j] && (C_sek[i][j] / C_glob[i][j] < 0.9 || C_sek[i][j] / C_glob[i][j] > 1.1))
                		{
                    			Err++;
                		}
            		}
        	}
        	printf("Błędy: %d (%.2f%%)\n", Err, (float)Err / (N * N) * 100);
        	
        	plik_out = fopen("Obl_równ.txt", "w");
        
        	if (plik_out == NULL)
        	{
            		printf("Blad otwarcia pliku \"wynik_row.txt\"\n");
        	}
        	else
       		{
            		for (int i = 0; i < N; i++)
            		{
                		for (int j = 0; j < N; j++)
                		{
                    			fprintf(plik_out, "%f ", C_glob[i][j]);
                		}
                		fprintf(plik_out, "\n");
            		}

            		fclose(plik_out);

        	}

        	plik_out = fopen("Obl_sekw.txt", "w"); 
        	if (plik_out == NULL)
        	{
            		printf("Blad otwarcia pliku \"wynik_sek.txt\"\n");
        	}
        	else
        	{
            		for (int i = 0; i < N; i++)
            		{
                		for (int j = 0; j < N; j++)
                		{
                    			fprintf(plik_out, "%f ", C_sek[i][j]);
                		}
                		fprintf(plik_out, "\n");
            		}

            		fclose(plik_out);
        	}

	}



	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
