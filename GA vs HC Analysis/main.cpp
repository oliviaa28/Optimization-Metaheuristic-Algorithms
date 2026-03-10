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

double decodificare_biti(vector<bool> biti, int start, double a, double b) {
    int decimal = 0;
    for (int i = 0; i < n; i++) {
        decimal = decimal * 2 + biti[start + i];
    }
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
    for (auto& indiv : populatie)
        indiv.cost = evaluare_solutie(functie, indiv.cromozom, a, b);
}

// roulette selection dimensiunea 5 si 10
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

// Tournament Selection pentru dimensiunea 30
Individ selectie_tournament(const vector<Individ>& pop, int tournament_size = 3) {
    uniform_int_distribution<int> dist(0, pop.size() - 1);
    
    Individ best = pop[dist(MT1)];
    for (int i = 1; i < tournament_size; i++) {
        Individ candidate = pop[dist(MT1)];
        if (candidate.cost< best.cost) best = candidate;
    }
    return best;
}

vector<Individ> selectie_tournament_pool(const vector<Individ>& pop, int new_size, int tournament_size = 3) {
    vector<Individ> selected;
    selected.reserve(new_size);
    
    for (int i=0; i<new_size; i++) 
        selected.push_back( selectie_tournament(pop,tournament_size) );
    
    return selected;
}

// Crossover cu 2 puncte 
void crossover(vector<Individ>& pop, double Pc) {
    if (pop.size() < 2) return;
    uniform_real_distribution<double> probabilitate(0.0, 1.0);
    uniform_int_distribution<int> punct(0, L - 1);

    // Amestecam perechile pentru diversitate
    vector<int> indices(pop.size());
    for (int i = 0; i < (int)pop.size(); i++) indices[i] = i;
    shuffle(indices.begin(), indices.end(), MT1);

    for (int i = 0; i + 1 < (int)pop.size(); i += 2) {
        if (probabilitate(MT1) < Pc) {
            int idx1 = indices[i],idx2 = indices[i+1];
            int p1 = punct(MT1), p2 =punct(MT1);

            if (p1 > p2) swap(p1, p2);
            if ( p1==p2 && p2< L-1 ) p2++;
            
            for ( int j=p1; j<=p2; j++ )  swap( pop[idx1].cromozom[j] , pop[idx2].cromozom[j] );
        }
    }
}

// Crossover Uniform pt Rastrigin D=30
void crossover_uniform(vector<Individ>& pop, double Pc){
    if (pop.size() <2) return;
    uniform_real_distribution<double> prob(0.0, 1.0);
    
    vector<int> indices( pop.size() );
    for (int i=0; i< (int)pop.size(); i++) indices[i] = i;
    shuffle( indices.begin(), indices.end(), MT1);

    for (int i= 0; i+1 <(int)pop.size() ; i+=2 ){
        if (prob(MT1) < Pc) {
            int idx1 = indices[i], idx2 = indices[i + 1];
            for (int j = 0; j < L; j++) {
                if (prob(MT1) < 0.5)             // Fiecare bit: 50% pr de swap
                    swap(pop[idx1].cromozom[j], pop[idx2].cromozom[j]);  
            }
        }
    }
}

void mutatie(vector<Individ>& pop, double Pm) {
    uniform_real_distribution<double> probabilitate(0.0, 1.0);
    double Pm_efectiv = max(Pm, 1.0 / L); //in functie de lungimea cromozomului
    
    for (int i = 0; i < (int)pop.size(); i++) {
        for (int j = 0; j < L; j++) {
            if (probabilitate(MT1) <Pm_efectiv)
                pop[i].cromozom[j]= !pop[i].cromozom[j];}
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

int calculeaza_elitism(int functie, int D, int marime_populatie) {
    int nr_elita = 0;
    if (D <= 10) {
        if (functie == 1) {
            nr_elita=(int)(0.15*marime_populatie);  // 15% elita
            nr_elita=max (nr_elita, 5);
        }
        else {
            nr_elita=(int)(0.12 * marime_populatie);  // 12% elita
            nr_elita= max(nr_elita, 4);
        }
    }
    // D=30 elitism redus 
    else {
        if (functie == 1) {// De Jong 10%
            nr_elita = (int)(0.10 *marime_populatie);
            nr_elita=max(nr_elita, 3);
        }else
         if (functie==3)  // RASTRIGIN 
             nr_elita =max(2, (int)(0.05* marime_populatie)) ; //  5%!!! 
        else {
            nr_elita =(int)(0.08* marime_populatie);
            nr_elita =max(nr_elita, 2);   }
    } 
    return nr_elita;
}

double algoritm_genetic_clasic(int functie, double a, double b, int marime_populatie,double prob_crossover, double prob_mutatie,int max_generatii)
{
    // Generare populatie initiala
    vector<Individ> populatie = genereaza_populatie(marime_populatie);
    evalueaza_populatie(populatie, functie, a, b);
    
    int nr_elita = calculeaza_elitism(functie, D, marime_populatie);
    vector<Individ> best_elita = gaseste_elita(populatie, nr_elita);
    double best_global_cost = best_elita[0].cost; 

    for (int gen = 0; gen < max_generatii; gen++) {
        
        int nr_copii=marime_populatie - nr_elita;
        vector<Individ> copii;

        if (D <= 10) 
            copii=selectie_roulette(populatie, nr_copii);
        else {
            int nr;
            if (functie == 3) nr=5;
            else nr=3;
            copii=selectie_tournament_pool(populatie, nr_copii, nr);
        }

      //  crossover(copii, prob_crossover);
      if (functie == 3 && D == 30)  crossover_uniform(copii, prob_crossover);  // Pentru Rastrigin D=30
      else                           crossover(copii, prob_crossover);          // Pentru restul

        mutatie(copii, prob_mutatie);
        evalueaza_populatie(copii, functie, a, b);
        
        //noua populație = elita +copii
        vector<Individ> populatie_noua;
        populatie_noua.insert( populatie_noua.end(), best_elita.begin(), best_elita.end() );
        populatie_noua.insert( populatie_noua.end(), copii.begin(), copii.end() );
        
        // Actualizam elita pt următoarea gen.
        best_elita = gaseste_elita(populatie_noua, nr_elita);
        
        if (best_elita[0].cost < best_global_cost)
            best_global_cost = best_elita[0].cost;
        
        populatie = populatie_noua;
    }
    return best_global_cost;
}

struct GA_param { //parametrii pt ga
    double Pc;          // Prob. crossover
    double Pm;          // Pr mutatie
    int pop_size;       // marime populatie 
    int max_gen;        // Nr maxim gen.
};

GA_param parametrii_dejong = {0.85, 0.012, 200, 2500};
GA_param parametrii_schwefel={0.85, 0.025, 200, 3000};
GA_param parametrii_rastrigin={0.80, 0.035, 200, 3000};
GA_param parametrii_michalewicz ={0.85, 0.045, 200, 3000};

void ruleaza_test_ag(int functie, double a, double b, const string& nume_functie, const string& nume_fisier) {
    GA_param parametrii;
    if (functie == 1)  parametrii=parametrii_dejong;
    else if (functie == 2) parametrii =parametrii_schwefel;
    else if (functie == 3) parametrii =parametrii_rastrigin;
    else if (functie == 4) parametrii =parametrii_michalewicz;

    ofstream fout(nume_fisier);
    vector<int> dimensiuni ={5,10, 30};
    int numar_rulari=40;
    
    for(int D_curent : dimensiuni ){
         D=D_curent;
         reprez_solutie(a, b);
        fout<<"\n________________________________________\n"<< "DIMENSIUNE D= " << D_curent << "\n"<< "__________________________________________\n\n";  
        cout<<" D="<<D_curent<< ":";

        if(functie ==3 && D_curent == 30){
          GA_param parametrii_rastrigin={0.75, 0.030, 350, 7000};
          parametrii =parametrii_rastrigin;
        }
        
       for (int rulare=1; rulare <= numar_rulari; rulare++) {

            auto start = chrono::steady_clock::now();
            double val_best=algoritm_genetic_clasic( functie, a, b, parametrii.pop_size, parametrii.Pc, parametrii.Pm, parametrii.max_gen);
            auto end = chrono::steady_clock::now();
            double timp_executie= chrono::duration<double> (end- start).count();
            
            fout<<"Rulare "<< rulare<< ": " << fixed << setprecision(5)<< val_best << "  (timp: " << fixed << setprecision(5) << timp_executie << "s)\n";
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
     executie_DeJong();
    executie_Schwefel();
    executie_Rastrigin();
     executie_Michalewicz();

    auto end_total=chrono::steady_clock::now();
    double timp_total=chrono::duration<double>(end_total - start_total).count();

    cout << "\n________________________________________________________________\n";
    cout << " Rulare completa\n" << "Timp total:" << fixed << setprecision(2) << timp_total <<" secunde (" << timp_total/60.0 << " minute)\n";
    return 0;
}