#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <random> 
#include <ctime>
#include <iomanip>
#include <chrono>
#include <thread>
#include <algorithm>
#include <numeric>
#include <omp.h> //pt paralelizare 
#include <limits> // Pentru numeric_limits

using namespace std; 
mt19937 MT0, MT1;
//#define MUT_RATE 0.2
#define CROSSOVER_RATE 0.4

// Permutare vector<int>
// Matrice vector<vector<int>>

struct dateInstanta {
    int n=0;            
    vector<vector<int>> flux;
    vector<vector<int>> distanta;
    int POP_SIZE = 200; 
}; 

void initializare_generator() {
    mt19937_64 helper;
    uint64_t base = time(nullptr) + clock() + hash<std::thread::id>{}(std::this_thread::get_id());
    helper.seed(base);
    helper.discard(10000); 
    
    MT0.seed(helper());
    MT1.seed(helper());
}

int nr_random(int min, int max){
    if(min>max) swap(min, max);
    uniform_int_distribution<> dist(min, max);
    return dist(MT1);
}

void citesteFisier(const string& fisier, dateInstanta& date) {
    ifstream fisier_in(fisier);
    int numar;
    fisier_in >> date.n;    //dimensiunea problemei

    date.flux.resize(date.n);//n *n
    
    for (int i = 0; i < date.n; i++) {
        date.flux[i].clear(); // goleste randul a anterior, daca exista 
        for (int j = 0; j < date.n; j++) {
            fisier_in>> numar;
            date.flux[i].push_back(numar);
        }
    }

    date.distanta.resize(date.n);
    
    for (int i=0; i <date.n; i++) 
    {
        date.distanta[i].clear();
        for (int j=0; j< date.n; j++) {
            fisier_in >>numar;
            date.distanta[i].push_back(numar);}
    }

    fisier_in.close();
}

//___________________________ QAP___________________________________________

 //fitess: cost mai mic = Fitness mai bun
double fitness(const vector<int>& p, const dateInstanta& date) {
    double sum= 0.0;
    for(int i=0; i<date.n; i++)
        for(int j=0; j<date.n ; j++)
           // F[i][j] * D[p[i]][p[j]]
            sum+= (double)date.flux[i][j] *(double) date.distanta[p[i]][p[j]];
    
    return sum;
}
 
//____________________________________________________ AG __________________________________________

vector<double> calc_costuri_pop( const vector<vector<int>>& populatie, const dateInstanta& data){
    int popSize=(int)populatie.size();
    vector<double> cost(popSize);

    #pragma omp parallel for
    for (int i=0; i<popSize; i++)
        cost[i]= fitness(populatie[i], data);
    
    return cost;
}

vector<vector<int>> initializeazaPop(const dateInstanta& data) {
    vector<vector<int>> populatie;

    for (int i= 0; i < data.POP_SIZE; i++) {
        vector<int> individ(data.n);

        iota(individ.begin(), individ.end(), 0);

        shuffle(individ.begin(), individ.end(), MT1);
        populatie.push_back(individ);
    }
    return populatie;
}


vector<vector<int>> selectieTurneu(const vector<vector<int>>& populatie, const vector<double>& cost, const dateInstanta& data) {
  
    int popSize=populatie.size(), k_turneu= 10; 
      vector<vector<int>> parinti(popSize);

    #pragma omp parallel
    {

        unsigned seed= (unsigned) time(nullptr)+(unsigned)clock() +(unsigned)(omp_get_thread_num() *1234567);
        mt19937 rng(seed);
        uniform_int_distribution<int> x(0, popSize -1);

    #pragma omp for  
    for (int i=0; i< popSize; i++) { //pt ca foloim threaduri, generam cate un seed diferit pt fiecaare

      int best_idx= x(rng);
      double best_fit= cost[best_idx];
       
        for (int j=1; j < k_turneu; j++) {
            int idx= x(rng);
            double fit= cost[idx];
            
            if (fit < best_fit) { 
                best_fit= fit;
                best_idx= idx;
            }
        }
      //  parinti.push_back( populatie[best_idx]);
      parinti[i]= populatie[best_idx];
    }
}
    return parinti;
}

vector<int> crossover(const vector<int>& p1, const vector<int>& p2, const dateInstanta& data, mt19937& rng) {
    uniform_real_distribution<double> prob(0.0, 1.0);
    uniform_int_distribution<int> dist_pos(0, data.n-1);

    if (prob(rng) >CROSSOVER_RATE)
        return p1;

    vector<int> copil(data.n, -1);
    vector<bool> vizitat(data.n, false);

    int t1= dist_pos(rng), t2= dist_pos(rng);
    if (t1 > t2) swap(t1, t2);

    for (int i=t1; i <= t2; i++){
        copil[i]= p1[i];
        vizitat[p1[i]]= true;
       }

    int c=0;
    for (int i=0; i < data.n; i++){
        if (copil[i] == -1){
            while (vizitat[ p2[c]] ) c++;
            copil[i] = p2[c++];}
     }

    return copil;
}

void mutatieSwap(vector<int>& cromozom, const dateInstanta& data, mt19937& rng, double mut_rate) {
    uniform_real_distribution<double> prob(0.0, 1.0);
    uniform_int_distribution<int> dist_pos(0, data.n - 1);

    if (prob(rng)< mut_rate){
        int a=dist_pos(rng), b= dist_pos(rng);
        while(a == b)
              b=dist_pos(rng);

        swap(cromozom[a], cromozom[b]);
    }
}

//________________________________________________________________________ Simulated Annealing (SA) ______________________________________________________________________________________________
double estimare_c0( vector<int> curent , double f_curent, const dateInstanta& date ){
// m1 = Let m1 be the total number of transitions proposed that improves strictly the value of objective function f<=0
// m2 =  let m2 be the number of other (indreasing) proposed transitions  f>0
// m0= m1+m2 be the total number of proposed transitions 
//  f=  media inratautirilor

    double c0 = 0.1;
	double suma = 0.0; //sum apentru f pozitive
	double acceptance_rate = 0.80 ;

	int m0, m1=0, m2=0 ;
	m0= 200; // numarul de tranzitii propuse

	uniform_int_distribution<int> dist(0, date.n-1);

    for( int i=0; i<m0; ++i){
        //aici generam vecinul prin swap
        vector<int> vecin= curent;
        int poz1, poz2;
        poz1 = dist(MT1);

        do{ poz2= dist(MT1); } while(poz1 == poz2);

        swap(vecin[poz1], vecin[poz2]);

        double f_vecin= fitness(vecin, date);
        if (f_vecin <= f_curent)
           m1++;
        else {
            m2++;
            suma+=(f_vecin - f_curent);
        }
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

vector<int> mutatie(const vector<int> &permutare, const dateInstanta& date){
    vector<int> vecin = permutare;

    if(date.n <2 )return vecin;
	uniform_int_distribution<int> dist(0, date.n - 1);

    int poz1= dist(MT1), poz2;
    do{ poz2= dist(MT1); } while(poz1 == poz2);

     swap(vecin[poz1], vecin[poz2]);
	return vecin;
}
//pt a face timpul de rulare mai scurt
double modificare_cost_schimb( const vector<int>& sol, int i, int j,const dateInstanta& date){
    double schimbare= 0.0;

    for (int k=0; k < date.n; k++) {
        if (k==i || k==j) continue;

        schimbare+=date.flux[i][k] *(date.distanta[sol[j]][sol[k]] -date.distanta[sol[i]][sol[k]]);
        schimbare+=date.flux[j][k] *(date.distanta[sol[i]][sol[k]] -date.distanta[sol[j]][sol[k]]);

        schimbare+=date.flux[k][i] *(date.distanta[sol[k]][sol[j]] -date.distanta[ sol[k]][sol[i]]);
        schimbare+=date.flux[k][j] *(date.distanta[sol[k]][sol[i]] -date.distanta[sol[k]][sol[j]]);
    }

    return schimbare;
}

vector<int> SA(vector<int> permutare, const dateInstanta& date )
{
    double alpha=0.99, Tmin=1e-6;

    int Lk= 200+10*date.n;

    vector<int> curent= permutare;
    double f_curent= fitness(curent, date);

    vector<int> best=curent;
    double f_best=f_curent;

    double ck= estimare_c0(curent, f_curent, date);
    int fara_imb= 0;

    uniform_real_distribution<double> prob(0.0, 1.0);

    do {
        bool imbunatatire=false;

        for (int l=0; l< Lk; l++) {

            int i=nr_random(0, date.n-1), j;
            do { j= nr_random(0, date.n - 1); } while (i == j);

            double delta=modificare_cost_schimb(curent, i, j, date);

            swap(curent[i], curent[j]);

            // Metropolis 
            if (delta <= 0.0) 
                f_curent += delta;
             else if (ck >0.0 &&  prob(MT1) < exp(-delta / ck)) 
                f_curent +=delta;
            else // revenim la solutia veche
                swap(curent[i], curent[j]);
            
            if (f_curent < f_best) {
                f_best = f_curent;
                best=curent;
                imbunatatire= true;
            }
        }

        ck *= alpha;
        if (imbunatatire) fara_imb = 0;
        else fara_imb++;

    } while (ck > Tmin && fara_imb < 70);

    return best;
}

double GA (const dateInstanta& data, vector<int>& best_perm_global_out){
    double mut_rate= 0.15;
    int stagnari= 0;
    
    vector<vector<int>> populatie= initializeazaPop(data);
    double best_global_fitness= numeric_limits<double>::max();
    vector<int> celMaiBunGlobal;

    int nrGeneratii= 2000; //!!!!!!

   for (int gen=0; gen < nrGeneratii; gen++) 
{
    vector<double> cost_pop=calc_costuri_pop(populatie, data);

    int index_best_local=0;
    double best_local_val=cost_pop[0];

    bool imb_global=false; 

    for (int i=0; i <(int)populatie.size(); i++) {
        double cost=cost_pop[i];

        if (cost < best_local_val) {
            best_local_val= cost;
            index_best_local= i;
        }

        if (cost < best_global_fitness) {
            best_global_fitness= cost;
            celMaiBunGlobal= populatie[i];
            imb_global=true;
        }
    }

    if (imb_global) stagnari=0;
    else stagnari++;

     //   vector<double> cost_pop= calc_costuri_pop(populatie, data);
       vector<vector<int>> parinti= selectieTurneu(populatie, cost_pop, data);

        vector<vector<int>> nouaGeneratie;
        nouaGeneratie.resize(data.POP_SIZE); //alocam direct 
        nouaGeneratie[0]=populatie[index_best_local]; // Elitism

       #pragma omp parallel
      {
    unsigned seed= (unsigned)time(nullptr)+ (unsigned)clock()+ (unsigned)(omp_get_thread_num()*1234567);

    mt19937 rng(seed);
    uniform_int_distribution<int> dist_p(0, parinti.size() - 1);
    uniform_real_distribution<double> prob(0.0, 1.0);

    #pragma omp for
    for (int i=1; i < data.POP_SIZE; i++) {
         int p1= dist_p(rng),p2=dist_p(rng);
        vector<int> copil= crossover(parinti[p1], parinti[p2], data, rng);
        mutatieSwap(copil, data, rng, mut_rate);

        nouaGeneratie[i] = copil;   
        }
    }

       vector<double> cost_nou= calc_costuri_pop(nouaGeneratie, data);

     int index_worst=0;
    double worst_val=cost_nou[0];

    for (int i=1; i < data.POP_SIZE; i++){
       if (cost_nou[i] > worst_val) {
            worst_val = cost_nou[i];
             index_worst=i;
             }
        }

        int best1=0, best2=1, best3=2;
        auto better= [&](int a, int b){ return cost_pop[a] < cost_pop[b]; };

        if (!better(best1,best2)) swap(best1, best2);
        if (!better(best1,best3)) swap(best1, best3);
        if (!better(best2, best3)) swap(best2, best3);

        vector<int> best_sa=populatie[best1];
        double best_sa_cost=cost_pop[best1];

         vector<int> sol1=SA(populatie[best1], data);
         vector<int> sol3=SA(populatie[best3], data);
         vector<int> sol2=SA(populatie[best2], data);

        double c1= fitness(sol1, data);
        if (c1 < best_sa_cost) 
        { best_sa_cost = c1;
             best_sa = sol1; }

        double c2= fitness(sol2, data);
        if (c2 < best_sa_cost)
         { best_sa_cost= c2; 
            best_sa= sol2; }

        double c3= fitness(sol3, data);
        if (c3 < best_sa_cost) 
        { best_sa_cost=c3; 
          best_sa = sol3; }

        vector<int> rafinata_SA=best_sa;
        double cost_rafinata=best_sa_cost;

        if (cost_rafinata < worst_val){
            nouaGeneratie[index_worst]=rafinata_SA;

            if (cost_rafinata < best_global_fitness){
                best_global_fitness=cost_rafinata;
                celMaiBunGlobal= rafinata_SA;
            }
        }

        populatie=nouaGeneratie;
    }

    best_perm_global_out= celMaiBunGlobal;
    return best_global_fitness;
}

int main() {
    initializare_generator();
  cout << fixed << setprecision(5);

vector<string> fisiere = {//"instante_test/had14.txt", "instante_test/esc16a.txt", "instante_test/chr12a.txt",
    
   "instante_test/nug30.txt",
 //  "instante_test/lipa30b.txt",
  //  "instante_test/ste36a.txt",
//    "instante_test/tho30.txt",
  //  "instante_test/esc64a.txt", //ok
    
    "instante_test/sko81.txt","instante_test/tai100a.txt"
    //, "instante_test/tai256c.txt"
    };

vector<string> nume = { //"had14.txt", "esc16a.txt", "chr12a.txt",
   "nug30.txt ",
   //"ste36a.txt",    "lipa30b.txt", "tho30.txt",  "esc64a.txt",
     
     "sko81.txt", "tai100a.txt"
     //, "tai256c.txt"
     };


    for (int f=0; f < (int)fisiere.size(); f++) {

        dateInstanta d;
        citesteFisier(fisiere[f], d);

        string outname= "rezultate_test/rezultate_" + nume[f];
        ofstream out(outname);

        out <<"Instanta: "<< fisiere[f] << "\n n=" << d.n << " POP_SIZE=" << d.POP_SIZE<< "\n\n";
        cout <<"\n_______________" << fisiere[f]<< " ( n=" << d.n << ") ___________________\n";


        for (int r= 1; r <= 30; r++) {
            vector<int> best_perm;

            auto t_start= chrono::high_resolution_clock::now();
            double best_cost= GA(d, best_perm);
            auto t_end= chrono::high_resolution_clock::now();

            double timp =chrono::duration<double>(t_end - t_start).count();
            double cost_real= fitness(best_perm, d);

            out <<"Run "<< r<< "  cost=" <<fixed << setprecision(5)<< cost_real<< " timp=" << timp << "\n";
            cout << "Run " << r << ": cost=" << cost_real<< " timp=" << timp<< "\n";
        }

        out.close();
    }

    return 0;
}