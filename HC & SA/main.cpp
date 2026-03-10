/*
N= numarul total de vaori
n= lungimea= numarul de biti necesar pt areprezenta N valori
L= lg totala a unei solutii complete
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <random> 
#include <ctime>
#include <iomanip>
#include<chrono>

#define MAX 1000
#define MAX_SA 10000

using namespace std;

mt19937 MT0, MT1; // MT0 -generator initial

int D;             //nr de dimensiuni = [ 5, 10, 30 ]
int functie;	  // tipul functiei 
int d = 5;         //precizia
int L;
int N, n;

void reprez_solutie(double a, double b) {
	N = (b - a) * pow(10, d);
	n = ceil(log2(N));
	L = n * D;
}

void initializare_generator() {

	 MT0.seed(time(nullptr)); // pt rularea finala
	//MT0.seed(28);				//pt teste
	for (int r = 0; r < 1000; r++) { MT0(); }

	MT1.seed(MT0());
}

vector<bool> generare_sol_initiala() {

	vector<bool> biti(L);
	uniform_real_distribution<float> generator(0.0, 1.0); //valori reale uniform distribuite intre 0.0 si 1.0

	for (int i = 0; i < L; i++) {
		float valoare = generator(MT1);

		if (valoare < 0.5)
			biti[i] = 0;
		else
			biti[i] = 1;
	}

	return biti;
}
//____________________________________________________________________________________________________________________________________
double decodificare_biti(vector<bool> biti, int start, double a, double b) {
	int decimal = 0;
	double x_real;

	for (int i = 0; i < n; i++) {
		decimal = decimal * 2 + biti[start + i];
	}

	x_real = a + decimal * (b - a) / (pow(2, n) - 1);  //X_real = a + decimal(xbiti) * (b - a) / (2^n - 1)
	return x_real;
}

vector<double> decodificare_sir_biti(vector<bool>& biti, double a, double b) { //transformarea btilor in valori reale
	vector<double> biti_decodati(D);
	int start;
	for (int i = 0; i < D; i++) {
		start = i * n;
		biti_decodati[i] = decodificare_biti(biti, start, a, b);
	}
	return biti_decodati;
}

//__________________________________________________________________ Functii ______________________________________________________________________________________________

//functie==1
double DeJong(int n, vector<double> x) { // f1(x)=sum(x(i)^2), i=1:n, -5.12<=x(i)<=5.12. !!indicii corespund formulei matematice 
	double suma = 0;
	for (int i = 0; i < n; i++)
		suma += x[i] * x[i];

	return suma;
}

//functie==2
double Schwefel(int n, vector<double> x) {// f7(x)=sum(-x(i)·sin(sqrt(abs(x(i))))), i=1:n; -500<=x(i)<=500.

	double suma = 0;
	for (int i = 0; i < n; i++)
		suma += x[i] * sin(sqrt(fabs(x[i])));

	return  - suma; 
}

//functie==3
double Rastrigin(int n, vector<double> x) { //f6(x)=10·n+sum(x(i)^2-10·cos(2·pi·x(i))), i=1:n; -5.12<=x(i)<=5.12. 

	double suma = 0;
	for (int i = 0; i < n; i++)
		suma += x[i] * x[i] - 10 * cos(2 * (atan(1) * 4) * x[i]);	//atan(1)=a rctg(1)= pi/4

	return 10 * n + suma;
}

//functie==4
double Michalewicz(int n, vector<double> x) {// f12(x)=-sum( sin(x(i)) ·(sin( i· x(i)^2/pi ))^(2·m)), i=1:n, m=10  0<=x(i)<=pi.

	double suma = 0;
	int m = 10;
	for (int i = 1; i <= n; i++)
		suma += sin(x[i-1]) * pow(sin((i * x[i-1] * x[i-1]) / M_PI), 2 * m);

	return - suma;
}


double evaluare_solutie(int functie, vector<bool> biti , double a, double b) { //functia fitness !!!
	vector<double> biti_decodati = decodificare_sir_biti(biti, a, b);

	     if (functie == 1)  return DeJong(D, biti_decodati);
	else if (functie == 2)  return Schwefel(D, biti_decodati);
	else if (functie == 3)  return Rastrigin(D, biti_decodati);
	else if (functie == 4)  return Michalewicz(D, biti_decodati);
	else return 0;

}

//________________________________________________________________________ Hill Climbing  (HC) _______________________________________________________________________________________
double HC_Best_Impr(int functie, double a, double b)
{
    int iteratie = 0;
    int stagnari = 0;
    const int MAX_STAGNARI = 50;

    vector<bool> sol_best = generare_sol_initiala();
    double val_best = evaluare_solutie(functie, sol_best, a, b);

    do {
        bool local = true;
        vector<bool> sol_curenta = sol_best;
        double val_curenta = val_best;

        do {
            local = true;
            vector<bool> vecin_best = sol_curenta;
            double cea_mai_buna = val_curenta;

            // caută cel mai bun vecin
            for (int i = 0; i < sol_curenta.size(); i++) {
                sol_curenta[i] = !sol_curenta[i];
                double val_vecin = evaluare_solutie(functie, sol_curenta, a, b);

                if (val_vecin < cea_mai_buna) {
                    vecin_best = sol_curenta;
                    cea_mai_buna = val_vecin;
                    local = false;
                }
                sol_curenta[i] = !sol_curenta[i];
            }

            if (!local) {
                sol_curenta = vecin_best;
                val_curenta = cea_mai_buna;
            }

        } while (local == false);

        if (val_curenta < val_best) {
            sol_best = sol_curenta;
            val_best = val_curenta;
            stagnari = 0;
        } 
        else stagnari++;

        if (stagnari > MAX_STAGNARI || (rand() % 100 )< 5) {
            sol_curenta = generare_sol_initiala();
            val_curenta = evaluare_solutie(functie, sol_curenta, a, b);

            if (val_curenta < val_best) {
                sol_best = sol_curenta;
                val_best = val_curenta;
            }
            stagnari = 0;
        }
        iteratie++;
    } while (iteratie < MAX);
    return val_best;
}

double HC_First_Impr(int functie, double a, double b)
{
    int iteratie = 0, stagnari = 0;
    const int MAX_STAGNARI = 50;

    vector<bool> sol_best = generare_sol_initiala();
    double val_best = evaluare_solutie(functie, sol_best, a, b);

    do {
        bool local = true;
        vector<bool> sol_curenta = sol_best;
        double val_curenta = val_best;

        do {
            local = true;
            for (int i = 0; i < sol_curenta.size(); i++) {
                sol_curenta[i] = !sol_curenta[i];
                double val_vecin = evaluare_solutie(functie, sol_curenta, a, b);

                // opreste laprima imbunatatire gasita 
                if (val_vecin < val_curenta) {
                    val_curenta = val_vecin;
                    local = false;
                    break;
                }
                sol_curenta[i] = !sol_curenta[i];
            }
        } while (local == false);

        // actualizare globala 
        if (val_curenta < val_best) {
            val_best = val_curenta;
            sol_best = sol_curenta;
            stagnari = 0;
        } else stagnari++;

      
        if (stagnari > MAX_STAGNARI || (rand() % 100 )< 5) 
        {
            sol_curenta = generare_sol_initiala();
            val_curenta = evaluare_solutie(functie, sol_curenta, a, b);
            if (val_curenta < val_best)
             {
                sol_best = sol_curenta;
                val_best = val_curenta;
            }
            stagnari = 0;
        }

        iteratie++;
    } while (iteratie < MAX);
    return val_best;
}

double HC_Worst_Impr(int functie, double a, double b)
{
    int iteratie = 0, stagnari = 0;
    const int MAX_STAGNARI = 50;

    vector<bool> sol_best = generare_sol_initiala();
    double val_best = evaluare_solutie(functie, sol_best, a, b);

    do {
        bool local = false;
        vector<bool> sol_curenta = sol_best;
        double val_curenta = val_best;

        do {
            local = true;
            vector<bool> vecin_best = sol_curenta;
            double cea_mai_rea = -1e9;
            bool gasit = false;

            for (int i = 0; i < sol_curenta.size(); i++) {
                sol_curenta[i] = !sol_curenta[i];
                double val_vecin = evaluare_solutie(functie, sol_curenta, a, b);

                if (val_vecin < val_curenta && (!gasit || val_vecin > cea_mai_rea)) {
                    vecin_best = sol_curenta;
                    cea_mai_rea = val_vecin;
                    gasit = true;
                    local = false;
                }
                sol_curenta[i] = !sol_curenta[i];
            }

            if (!local) {
                sol_curenta = vecin_best;
                val_curenta = cea_mai_rea;
            }

        } while (local == false);

        //actiualizam solutia globala
        if (val_curenta < val_best) {
            sol_best = sol_curenta;
            val_best = val_curenta;
            stagnari = 0;
        } else stagnari++;

        // restart dacă algoritmul stagnează sau cu o mică probabilitate
        if (stagnari > MAX_STAGNARI || rand() % 100 < 5) {
            sol_curenta = generare_sol_initiala();
            val_curenta = evaluare_solutie(functie, sol_curenta, a, b);
            if (val_curenta < val_best) {
                sol_best = sol_curenta;
                val_best = val_curenta;
            }
         stagnari = 0;
        }

    iteratie++;

    } while (iteratie < MAX);

    return val_best;
}

//________________________________________________________________________ Simulated Annealing (SA) ______________________________________________________________________________________________
double estimare_c0( vector<bool> curent , double f_curent , int functie, double a, double b){
// m1 = Let m1 be the total number of transitions proposed that improves strictly the value of objective function f<=0
// m2 =  let m2 be the number of other (indreasing) proposed transitions  f>0
// m0= m1+m2 be the total number of proposed transitions 

//  f=  media inratautirilor

    double c0 = 0.1;
	double suma = 0.0; //sum apentru f pozitive
	double acceptance_rate = 0.80 ;

	int m0, m1 = 0, m2 = 0 ;
	m0= 200; // numarul de tranzitii propuse

	uniform_int_distribution<int> dist(0, L-1);

	while ( m0){

		int i = dist (MT1); 
		vector<bool> vecin = curent;
		vecin[i] = !vecin[i];
		double f_vecin = evaluare_solutie(functie, vecin, a, b);

		  if ( f_vecin - f_curent <= 0)
		      m1++;
		else {
			m2++;
			suma += (f_vecin - f_curent);}

		m0--;
	}

	if (m2 == 0) return 1e-3; //evitam impartirea la 0

	double media_inrautatirilor , log_arg; 
	//c0 = media_inrautatirilor / ( log(m2 /( m2* acceptance_rate - m1*( 1 - acceptance_rate) ) ) );
	//acceptance_rate = ( m1+ m2 * exp( -(media_inrautatirilor / c0 )) ) /( m1+m2 );

    media_inrautatirilor = suma / m2; 

	int deimpartit = m2* acceptance_rate - m1*( 1.0 - acceptance_rate);
	if ( deimpartit <=0) deimpartit = 1e-3;

	log_arg = (double)m2 / (double)deimpartit;
	if ( log_arg <= 1) log_arg = 1 + 1e-3; 
     
	c0 = media_inrautatirilor / ( log(log_arg) );

	return c0;
}

vector<bool> mutatie(vector<bool> curent , int functie){
	uniform_int_distribution<int> dist(0, L - 1);

	
	if (functie == 1){
		int i= dist(MT1);
		curent[i] = !curent[i];
	}
	else if (functie ==2 ){
		int k= 1+ (dist(MT1)%7); //modificam intre 1 si 3 biti
		for ( int j=0; j< k; j++){
			int i= dist(MT1);
			curent[i] = !curent[i];
		}
	}
	else if (functie ==3 ){
		int k= 1+ (dist(MT1)%5); //  1 si 5 biti
		for ( int j=0; j< k; j++){
			int i= dist(MT1);
			curent[i] = !curent[i];
		}
	}
	else {
		int k= 1+ ( dist(MT1) %7 ); // 1 si 7 biti
		for ( int j=0; j< k; j++){
			int i= dist(MT1);
			curent[i] = !curent[i];
		}
	}

	return curent;
}

double SA(int functie, double a, double b){

	//initializare 
	//Lk = the number of transitions generated by some iteration k . lungimea lantului la fiecare temp
	// ck=the value of the temperature paraemter
	double alpha = 0.98 ;
	double Tmin = 1e-6;

	int i, k, max_iteratii= MAX_SA;
	int Lk  =  200 + 10 * D; 

	vector<bool> curent = generare_sol_initiala();
	double f_curent = evaluare_solutie(functie, curent , a , b);

	vector<bool> best = curent;
	double f_best = f_curent; 

	double c0 = estimare_c0( curent , f_curent , functie, a, b);
	if (functie == 2) c0=1.0;
//	double c0 = 1e-3;// (0.001)
	double ck = c0;

	int fara_imbunatariri=0;
	bool imbunatatire= 0;

	uniform_real_distribution<double> prob(0.0, 1.0);

	do{
		bool imbunatatire= 0;
		for(int l = 0; l < Lk; l++) {
			//rand neighbours
			vector<bool> vecin = mutatie( curent , functie);
			double f_vecin = evaluare_solutie(functie, vecin, a, b);

			//Metropolis acceptance criterion
			if(f_vecin <= f_curent) { //pr =1
			   f_curent = f_vecin;
			    curent = vecin;
			} else {
                if ( ck > 0)
			       if( ( prob(MT1) < exp(- (f_vecin - f_curent) / ck) )) {
				  		 f_curent = f_vecin;
				    	 curent = vecin; }
			}

			if(f_curent < f_best) {
			    f_best = f_curent;
			    best = curent;
				imbunatatire =1;
			}	
		} 

		max_iteratii--;

		//geometric cooling 
		 ck *= alpha;
		 if (imbunatatire == 0) fara_imbunatariri++;
		 else fara_imbunatariri=0;

	} while(ck > Tmin && fara_imbunatariri < 66 && max_iteratii > 0);

	 return f_best;
}
//____________________________________________________________________________________________________________________________________________________________________________________________________

void executie_DeJong() {
    double a = -5.12, b = 5.12;
    int functie = 1;

 //   cout << "______________________________ Test Hill Climbing - De Jong ___________________________________\n";

    vector<int> dimensiuni = {5, 10, 30};

    for (int D_curent : dimensiuni) {

        D = D_curent;
        reprez_solutie(a, b); 
        cout << "\n=== Dimensiune: " << D << " ===\n";

        auto start = chrono::steady_clock::now();
        double val_best = HC_Best_Impr(functie, a, b);
        auto end = chrono::steady_clock::now();
        double timp_best = chrono::duration<double>(end - start).count();
        cout << fixed << setprecision(5);
        cout << "[BEST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_best << " in " << timp_best << " secunde.\n\n";

        auto start1 = chrono::steady_clock::now();
        double val_first = HC_First_Impr(functie, a, b);
        auto end1 = chrono::steady_clock::now();
        double timp_first = chrono::duration<double>(end1 - start1).count();
        cout << "[FIRST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_first << " in " << timp_first << " secunde.\n\n";

        auto start2 = chrono::steady_clock::now();
        double val_worst = HC_Worst_Impr(functie, a, b);
        auto end2 = chrono::steady_clock::now();
        double timp_worst = chrono::duration<double>(end2 - start2).count();
        cout << "[WORST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_worst << " in " << timp_worst << " secunde.\n\n";
    }

    cout << "_______________________________________________________________________________________________\n";

	// cout << "______________________________ Test SA - De Jong ___________________________________\n";

	/*
    auto t0 = chrono::steady_clock::now();
     double f_sa = SA(functie, a, b);
     auto t1 = chrono::steady_clock::now();
    double sec_sa = chrono::duration<double>(t1 - t0).count();
    cout << "[SA, dim 30] Valoare minima gasita: " << f_sa
         << "  | timp(min): " << sec_sa  << " s\n";

*/
}
void executie_Schwefel() {
    int functie = 2;
    double a = -500, b = 500;
    vector<int> dimensiuni = {5, 10, 30};

  //  cout << "______________________________ Test Hill Climbing - Schwefel ___________________________________\n";

    for (int D_curent : dimensiuni) {
        D = D_curent;
        d = 5;                 
        reprez_solutie(a, b);  

        cout << "\n=== Dimensiune: " << D << " ===\n";
/*
        auto start = chrono::steady_clock::now();
        double val_best = HC_Best_Impr(functie, a, b);
        auto end = chrono::steady_clock::now();
        double timp_best = chrono::duration<double>(end - start).count();
        cout << fixed << setprecision(d);
        cout << "[BEST Improvement, dim=" << D << "] Valoare minima gasita: " 
             << val_best << " in " << timp_best << " secunde. \n\n";

       auto start1 = chrono::steady_clock::now();
        double val_first = HC_First_Impr(functie, a, b);
        auto end1 = chrono::steady_clock::now();
        double timp_first = chrono::duration<double>(end1 - start1).count();
        cout << "[FIRST Improvement, dim=" << D << "] Valoare minima gasita: " 
             << val_first << " in " << timp_first << " secunde. \n\n";
*/
     auto start2 = chrono::steady_clock::now();
        double val_worst = HC_Worst_Impr(functie, a, b);
        auto end2 = chrono::steady_clock::now();
        double timp_worst = chrono::duration<double>(end2 - start2).count();
        cout << "[WORST Improvement, dim=" << D << "] Valoare minima gasita: " 
             << val_worst << " in " << timp_worst << " secunde. \n\n";
	
    }

    cout << "_________________________________________________________________________________________________\n";


	 /*
  
	 //   cout << "______________________________ Test SA - Schwefel ___________________________________\n";


     auto t0 = chrono::steady_clock::now();
     double f_sa = SA(functie, a, b);
     auto t1 = chrono::steady_clock::now();
     double sec_sa = chrono::duration<double>(t1 - t0).count();
     cout << "[SA , dim 5] Valoare minima gasita: " << f_sa
         << "  | timp: " << sec_sa  << " s\n";

*/
}
void executie_Rastrigin() {
    int functie = 3;
    double a = -5.12, b = 5.12;
    vector<int> dimensiuni = {5, 10, 30};

    cout << "______________________________ Test Hill Climbing - Rastrigin __________________________________\n";

    for (int D_curent : dimensiuni) {
        D = D_curent;
        d = 5;
        reprez_solutie(a, b);

        cout << "\n=== Dimensiune: " << D << " ===\n";

        auto start = chrono::steady_clock::now();
        double val_best = HC_Best_Impr(functie, a, b);
        auto end = chrono::steady_clock::now();
        double timp_best = chrono::duration<double>(end - start).count();
        cout << fixed << setprecision(d);
        cout << "[BEST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_best << " in " << timp_best << " secunde. \n\n";

        auto start1 = chrono::steady_clock::now();
        double val_first = HC_First_Impr(functie, a, b);
        auto end1 = chrono::steady_clock::now();
        double timp_first = chrono::duration<double>(end1 - start1).count();
        cout << "[FIRST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_first << " in " << timp_first << " secunde. \n\n";

        auto start2 = chrono::steady_clock::now();
        double val_worst = HC_Worst_Impr(functie, a, b);
        auto end2 = chrono::steady_clock::now();
        double timp_worst = chrono::duration<double>(end2 - start2).count();
        cout << "[WORST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_worst << " in " << timp_worst << " secunde. \n\n";
    }

    cout << "_________________________________________________________________________________________________\n";


	//  cout << "______________________________ Test SA - RASTRIGIN___________________________________\n";*/
/*
     auto t0 = chrono::steady_clock::now();
     double f_sa = SA(functie, a, b);
     auto t1 = chrono::steady_clock::now();
    double sec_sa = chrono::duration<double>(t1 - t0).count();
    cout << "[SA , dim 30 ] Valoare minima gasita: " << f_sa
         << "  | timp: " << sec_sa << " s\n"; */

}
void executie_Michalewicz() {
    int functie = 4;
    double a = 0.0, b = M_PI;
    vector<int> dimensiuni = {5, 10, 30};

    cout << "______________________________ Test Hill Climbing - Michalewicz ________________________________\n";

    for (int D_curent : dimensiuni) {
        D = D_curent;
        d = 5;
        reprez_solutie(a, b);

        cout << "\n=== Dimensiune: " << D << " ===\n";

        auto start = chrono::steady_clock::now();
        double val_best = HC_Best_Impr(functie, a, b);
        auto end = chrono::steady_clock::now();
        double timp_best = chrono::duration<double>(end - start).count();
        cout << fixed << setprecision(d);
        cout << "[BEST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_best << " in " << timp_best << " secunde. \n\n";

        auto start1 = chrono::steady_clock::now();
        double val_first = HC_First_Impr(functie, a, b);
        auto end1 = chrono::steady_clock::now();
        double timp_first = chrono::duration<double>(end1 - start1).count();
        cout << "[FIRST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_first << " in " << timp_first << " secunde. \n\n";

        auto start2 = chrono::steady_clock::now();
        double val_worst = HC_Worst_Impr(functie, a, b);
        auto end2 = chrono::steady_clock::now();
        double timp_worst = chrono::duration<double>(end2 - start2).count();
        cout << "[WORST Improvement dim=" << D << "] Valoare minima gasita: " 
             << val_worst << " in " << timp_worst << " secunde. \n\n";
    }

    cout << "_________________________________________________________________________________________________\n";

	// cout << "______________________________ Test SA - Michalewicz___________________________________\n";
/*
    auto t0 = chrono::steady_clock::now();
     double f_sa = SA(functie, a, b);
     auto t1 = chrono::steady_clock::now();
    double sec_sa = chrono::duration<double>(t1 - t0).count();
    cout << "[SA , dim 5] Valoare minima gasita: " << f_sa
         << "  | timp: " << sec_sa  << " s\n";   */

	
}

int main() {
    
    initializare_generator();
	cout << fixed << setprecision(5);

	int nr=30;
     while(nr){ executie_DeJong(); nr--;}

	  nr=30;
	//  while(nr){ executie_Schwefel(); nr--;}

	  nr=30;
	//  while(nr){ executie_Rastrigin(); nr--;}
	
	  nr=30;
	//  while(nr){ executie_Michalewicz(); nr--;}

	return 0;
}