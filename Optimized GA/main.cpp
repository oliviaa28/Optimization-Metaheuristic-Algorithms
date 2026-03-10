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

#include <omp.h>

using namespace std;

mt19937 MT0, MT1;
int D, functie, d=5, L, N, n;

void reprez_solutie(double a, double b) {
    N = (b - a) * pow(10, d);
    n = ceil(log2(N));
    L = n * D;
}

void initializare_generator() {
    mt19937_64 helper;
    uint64_t base = time(nullptr) + clock() + std::hash<std::thread::id>{}(std::this_thread::get_id());
    helper.seed(base);
    helper.discard(10000); //scapam de primele 10000
    
    MT0.seed(helper());
    MT1.seed(helper());
}

//____________________________________________________________ reprezentae solutii ____________________________________________________________________________________________________

vector<bool> binary_to_gray(vector<bool> binar){
    vector<bool> gray( binar.size()); 
    gray[0] = binar[0];    //primul bit ramane la fel

    for(int i=0; i <(int)binar.size()-1; i++)
        gray[i+1]= binar[i] ^ binar[i+1]; // XOR intre bitul curent si urmatorul
    
    return gray;
}

vector<bool> gray_to_binary(const vector<bool>& gray){
    vector<bool> binar(gray.size());
    binar[0]= gray[0];
    for( int i=1; i< (int)gray.size(); i++ )
   {  
       if(gray[i] == 0 )
         binar[i]= binar[i-1];
     else 
        binar[i]= !binar[i-1];
    }
    
        return binar;
}

double decodificare_biti(vector<bool> biti, int start, double a, double b) {
    vector<bool> segment(n);
    for (int i=0; i<n; i++)
        segment[i] = biti[start+i];
    
    segment= gray_to_binary(segment);

    int decimal = 0;
    for (int i = 0; i < n; i++) {
        decimal = decimal * 2 + segment[i];}

    double x_real = a + decimal * (b - a) / (pow(2, n) - 1); 
    return x_real;
}

vector<double> decodificare_sir_biti(vector<bool>& biti, double a, double b) {
    vector<double> biti_decodati(D);
    for (int i = 0; i < D; i++) {
        int start = i * n;
        biti_decodati[i] = decodificare_biti(biti, start, a, b);
    }
    return biti_decodati;
}

//________________________________________________________________________________________________________________________________________________________________
double DeJong(int n, vector<double> x) { 
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += x[i] * x[i];
    return suma;
}

double Schwefel(int n, vector<double> x) {
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += x[i] * sin(sqrt(fabs(x[i])));
    return -suma; 
}

double Rastrigin(int n, vector<double> x) { 
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]);
    return 10 * n + suma;
}

double Michalewicz(int n, vector<double> x) {
    double suma = 0;
    int m = 10;
    for (int i = 1; i <= n; i++)
        suma += sin(x[i-1]) * pow(sin((i * x[i-1] * x[i-1]) / M_PI), 2 * m);
    return -suma;
}

double evaluare_solutie(int functie, vector<bool> biti, double a, double b) {
    vector<double> biti_decodati = decodificare_sir_biti(biti, a, b);
    
    if (functie == 1)       return DeJong(D, biti_decodati);
    else if (functie == 2)  return Schwefel(D, biti_decodati);
    else if (functie == 3)  return Rastrigin(D, biti_decodati);
    else if (functie == 4)  return Michalewicz(D, biti_decodati);
    else return 0;
}

//________________________________________________________________________ HC ______________________________________________________________________________________________
double hill_climbing_best_improvement(int functie, double a, double b,
                                      vector<bool> solutie_ga, int max_iteratii = 100)
{
    vector<bool> sol_best = solutie_ga;
    double val_best = evaluare_solutie(functie, sol_best, a, b);

    int stagnari = 0;
    const int MAX_STAGNARI = 20;

    vector<bool> vecin_best;
    double val_vecin_best;

    for (int iter = 0; iter < max_iteratii && stagnari < MAX_STAGNARI; iter++) {

        bool gasit = false;
        vecin_best = sol_best;
        val_vecin_best = val_best;

#pragma omp parallel for
        for (int i = 0; i < L; i++) {
            vector<bool> copie = sol_best;
            copie[i] = !copie[i];

            double val = evaluare_solutie(functie, copie, a, b);

#pragma omp critical
            {
                if (val < val_vecin_best) {
                    val_vecin_best = val;
                    vecin_best = copie;
                    gasit = true;
                }
            }
        }

        if (gasit) {
            sol_best = vecin_best;
            val_best = val_vecin_best;
            stagnari = 0;
        } else stagnari++;
    }
    return val_best;
}

//__________________________________________________________________algoritm genetic ______________________________________________________________________________________________
struct Individ {
    vector<bool> cromozom;
    double cost;
};

vector<bool> genereaza_cromozom_random() {
    vector<bool> cromozom(L);
    uniform_int_distribution<int> dist_int(0, 1);
    for (int i = 0; i < L; i++) cromozom[i] = dist_int(MT1);
	
    return cromozom;
}

vector<Individ> genereaza_populatie(int pop_size) {
    vector<Individ> populatie(pop_size);
    for (int i = 0; i < pop_size; i++)   populatie[i].cromozom=genereaza_cromozom_random();
    
    return populatie;
}


void evalueaza_populatie(vector<Individ>& populatie, int functie, double a, double b) {
    #pragma omp parallel for
    for (int i = 0; i < (int)populatie.size(); i++)
        populatie[i].cost = evaluare_solutie(functie, populatie[i].cromozom, a, b);
}


vector<Individ> selectie_roulette(const vector<Individ>& pop, int new_size) {
    int P = (int)pop.size();
    vector<Individ> rez;
    rez.reserve(new_size);

    // Gasim min și max cost
    double min_cost = pop[0].cost, max_cost = pop[0].cost;
    for (const auto& ind : pop) {
        if (ind.cost < min_cost) min_cost = ind.cost;
        if (ind.cost > max_cost) max_cost = ind.cost;
    }

    double eps = 1e-12, range = max_cost - min_cost;

    // Daca toti au aproape acelasi cost, selectie uniforma
    if (range < eps) {
        uniform_int_distribution<int> rndIndex(0, P - 1);
        for (int k = 0; k < new_size; k++) {
            rez.push_back( pop[ rndIndex(MT1) ] );
        }
        return rez;
    }

    vector<double> scores(P), cumulat(P);
    double sum_scores = 0.0;

    for (int i = 0; i < P; i++) {
        scores[i] = (max_cost - pop[i].cost) + eps;
        sum_scores += scores[i];
    }

    // Construim functia cumulatulativa
    cumulat[0] = scores[0] / sum_scores;
    for (int i = 1; i < P; i++)
        cumulat[i]= cumulat[i-1] + scores[i]/sum_scores;
    cumulat[P-1] = 1.0;

    uniform_real_distribution<double> dist(0.0, 1.0);
    for (int k =0; k <new_size;k++) //selectia 
    {
        double r = dist(MT1);
        for (int j = 0; j < P; j++) {
            if (r <= cumulat[j]) {
                rez.push_back(pop[j]);
                break;
            }
        }
    }
    return rez;
}

//___________________________________________________________ Adaptive probabilities  ______________________________________________________________________________________________
struct StatFitness {
    double fmax;
    double favg;
};

StatFitness calc_stat( const vector<Individ>& pop){
    double fmax = pop[0].cost, suma = 0.0, favg;

    for (const auto& ind : pop) {
        if (ind.cost < fmax)  
            fmax=ind.cost;//la minimizare
        suma +=ind.cost;
    }
    favg= suma / pop.size();
    return {fmax, favg};
}

void crossover_adaptiv(vector<Individ>& pop, double k1, double k3, double fmax, double favg)
{       if(pop.size() < 2) return;
     uniform_real_distribution<double> rand01( 0.0, 1.0 );
    uniform_int_distribution<int> punct( 0, L - 1);

    vector<int> indices(pop.size());
    iota( indices.begin(), indices.end(),0 );
    shuffle( indices.begin(), indices.end(), MT1);

    double delta= fabs(fmax - favg);
    if (delta < 1e-12) delta = 1e-12; //sa evitam impartirea la 0

    for (int i=0; i+1<(int)pop.size(); i += 2){
        int a =indices[i], b =indices[i+1 ];
        double fprime =min(pop[a].cost, pop[b].cost); // la minimizare
        double pc;

        if (fprime <= favg)
            pc=k1 * fabs(fmax - fprime) / delta;
        else
            pc=k3;

        if (rand01(MT1) < pc) {
            int p1 = punct(MT1), p2 = punct(MT1);
            if (p1 > p2) swap(p1, p2);
            if (p1==p2 && p2<L-1) p2++;

            for (int j = p1; j <= p2; j++) swap(pop[a].cromozom[j], pop[b].cromozom[j]);
        }
    }
}
/*
void mutatie_adaptiva( vector<Individ>& pop, double k2 , double k4 ,double fmax,double favg)
{
    uniform_real_distribution<double> rand01(0.0, 1.0);
    double delta=fabs(fmax - favg);
    if (delta < 1e-12) delta = 1e-12; //sa evitam impartirea la 0

    for (auto& ind:pop) {//pentru fiecare individ
        double f = ind.cost, pm;

        if (f<=favg)
            pm=k2*fabs(fmax - f) / delta;
        else
            pm=k4;

        for (int j=0; j<L; j++)
            if (rand01(MT1) < pm)
                ind.cromozom[j]= !ind.cromozom[j];
    }
}
*/

void mutatie_ga(vector<Individ>& pop, double p_max) { //Rank Based Adaptive Mutation Probability
//p = p_max*(1- (r-1)/(population_size -1 ) ))
//p=mutation probability of a chromosone 
//p_max = maximum muattion probability 
//r=rank of a chromosone 
  /* sortam populatia dupa cost si dam rankuri :
  cel mai bun are rank= pop_size , cel mai rau rank 1
   apoi aplicam formula pentru fiecare individ  */

  int pop_size = (int)pop.size();
    vector<int> index(pop_size);
    double p;
    for (int i = 0; i < pop_size; i++) 
        index[i] = i;

    // sortam dupa cost
    sort(index.begin(), index.end(), [&](int a, int b){
        return pop[a].cost < pop[b].cost;
    });

    vector<int> rank(pop_size);
    for (int i = 0; i < pop_size; i++) {
        int idx= index[i];
        rank[idx]= pop_size - i;  
    }

    uniform_real_distribution<double> prob(0.0, 1.0);

    for (int i = 0; i < pop_size; i++) {
        p =p_max * (1.0 -(double)(rank[i] - 1) / (double)(pop_size - 1) );
        for (int j =0; j<L; j++) {
            if ( prob(MT1) < p)
                pop[i].cromozom[j] = !pop[i].cromozom[j];
        }
    }
}

// gasim cei mai buni k indivizi
vector<Individ> gaseste_elita(const vector<Individ>& pop, int k) {
    vector<Individ> copie = pop;
    sort( copie.begin(), copie.end(), [](const Individ& a, const Individ& b) {
        return a.cost < b.cost;
    });
    
    vector<Individ> elita;
    for ( int i = 0; i<k && i<(int)copie.size(); i++ )
	 elita.push_back(copie[i]);
    
    return elita;
}

int calc_elita(int pop_size, int functie) {
    if (D == 30) 
        if (functie== 2 || functie==3 || functie==4) 
            return max(2, pop_size / 20); // 5%
        
     return max(3, pop_size / 10); //10%
}
//____________________________________________________________ Algoritm genetic  ______________________________________________________________________________________________
double algoritm_genetic(int functie, double a, double b, int marime_populatie,double prob_crossover, double prob_mutatie,int max_generatii)
{
    // Generare populatie initiala
    vector<Individ> populatie = genereaza_populatie(marime_populatie);
    evalueaza_populatie(populatie, functie, a, b);
    
    int nr_elita =calc_elita (marime_populatie, functie);
    vector<Individ> best_elita = gaseste_elita(populatie, nr_elita);
    double best_global_cost = best_elita[0].cost; 

    auto stats=calc_stat(populatie);
    double fmax =stats.fmax,favg =stats.favg;

    for (int gen = 0; gen < max_generatii; gen++) {
        auto stats = calc_stat(populatie);
        int nr_copii= marime_populatie -nr_elita;
        vector<Individ> copii=selectie_roulette(populatie, nr_copii);

        crossover_adaptiv(copii,1.0, 0.5, fmax,favg);
         mutatie_ga(copii, prob_mutatie);
        evalueaza_populatie(copii, functie, a, b);

        vector<Individ> populatie_noua;
        populatie_noua.insert( populatie_noua.end(),best_elita.begin(), best_elita.end());
        populatie_noua.insert( populatie_noua.end(),copii.begin(), copii.end() );
        
        best_elita = gaseste_elita(populatie_noua, nr_elita);
        
        if (best_elita[0].cost < best_global_cost)
            best_global_cost = best_elita[0].cost;
        
        populatie = populatie_noua;
    }
   // return best_global_cost;
   return hill_climbing_best_improvement(functie, a, b, best_elita[0].cromozom);
}

//____________________________________________________ EXECUTIIII_______________________________________________________________________

struct GA_param { double Pc;  double Pm;   int pop_size;  int max_gen; };
GA_param universal = {0.85, 0.015, 100, 1000 };

void ruleaza_test_ag(int functie, double a, double b, const string& nume_functie, const string& nume_fisier) {
    GA_param parametrii;
    parametrii = universal;

    ofstream fout(nume_fisier);
 // vector<int> dimensiuni ={5, 10, 30};
    vector<int> dimensiuni ={30};
    int numar_rulari=30;
    
    for(int D_curent : dimensiuni ){
         D=D_curent;
         reprez_solutie(a, b);
        fout<<"\n________________________________________\n"<< "DIMENSIUNE D= " << D_curent << "\n"<< "__________________________________________\n\n";  
        cout<<" functie= "<<functie<<" | D = "<<D_curent<< ":";
        
        if (D_curent == 30) {
            if (functie == 2) { // Schwefel
                parametrii.pop_size = 300;
                parametrii.max_gen  = 3000;
                parametrii.Pm       = 0.03; // mai mare
            } else if (functie == 3) { // Rastrigin
                parametrii.pop_size = 300;
                parametrii.max_gen  = 3500;
                parametrii.Pm       = 0.025;
            } else { // DeJong / Michalewicz
                parametrii.pop_size = 200;
                parametrii.max_gen  = 2000;
                parametrii.Pm       = 0.02;
            }
        }
        
       for (int rulare=1; rulare <= numar_rulari; rulare++) {

            auto start = chrono::steady_clock::now();
            double val_best=algoritm_genetic( functie, a, b, parametrii.pop_size, parametrii.Pc, parametrii.Pm, parametrii.max_gen);
            auto end = chrono::steady_clock::now();
            double timp_executie= chrono::duration<double> (end- start).count();
            
            fout<<"Rulare "<< rulare<< ": " << fixed << setprecision(5)<< val_best << "  (timp: " << fixed << setprecision(5) << timp_executie << "s)\n";
            cout<< fixed << setprecision(5)<< val_best<<endl;
           if (rulare% 5==0) cout<<" | functie= "<<functie<<" | rulare " << rulare << "..."<< flush; //ca sa avem putin de feedback in terminal
       }
    }   
    fout.close();
}

void executie_DeJong() {
    ruleaza_test_ag(1, -5.12, 5.12, "De Jong", "rezultate_DeJong.txt");
}

void executie_Schwefel() {
    ruleaza_test_ag(2, -500.0, 500.0, "Schwefel", "rezultate_Schwefel.txt");
}

void executie_Rastrigin() {
    ruleaza_test_ag(3, -5.12, 5.12, "Rastrigin", "rezultate_Rastrigin.txt");
}

void executie_Michalewicz() {
    ruleaza_test_ag(4, 0.0, M_PI, "Michalewicz", "rezultate_Michalewicz.txt");
}

int main() {
    initializare_generator();
    cout << fixed << setprecision(5);
   auto start_total = chrono::steady_clock::now();

 // executie_DeJong();
  // executie_Schwefel();
   executie_Rastrigin();
  //executie_Michalewicz();

    auto end_total=chrono::steady_clock::now();
    double timp_total=chrono::duration<double>(end_total - start_total).count();

    cout<< "\n________________________________________________________________\n";
    cout<< " Rulare completa\n" << "Timp total:" << fixed << setprecision(2) << timp_total <<" secunde (" << timp_total/60.0 << " minute)\n";
    return 0;
}