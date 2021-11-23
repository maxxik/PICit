//-----------------------------------------------------------------//
//       PICit! : flexible 1d3v PIC/MCC simulation code            //
//-----------------------------------------------------------------//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <random>
#include <numeric>
#include <functional>
#include <algorithm>
#include <chrono>

#if defined(__INTEL_COMPILER)
# include "vectorclass/version1/vectorclass.h" 
#else
# include "vectorclass/version2/vectorclass.h" 
#endif

#define SQR(x) ((x) * (x))

using namespace::std;

// constants
#include "consts.h"

// simulation parameters
#include "params.h"

// collision cross sections

int              N_CS [N_SPECIES];                           // number of processes for each species
string           cs_names [N_SPECIES];                       // cs dataset names for documentation
typedef vector<float>  cross_section;                        // cross section vector
cross_section    sigma [N_SPECIES];                          // set of cross section arrays
cross_section    sigma_tot [N_SPECIES][N_TARGET];            // total macroscopic cross section of electrons
double           max_coll_probability [N_SPECIES];           // maximum collision probabilit for null-collisions
double           max_coll_freq0 [N_SPECIES][N_TARGET];       // maximum collision frequency for unit density
vector<int>      sigma_type [N_SPECIES];                     // collision type identifiers
vector<int>      sigma_target [N_SPECIES];                   // target particle identifiers
vector<double>   sigma_threshold [N_SPECIES];                // collision type energy thresholds
vector<int>      sigma_Nproduct [N_SPECIES];                 // number of collision products 
vector<int>      sigma_product [N_SPECIES][4];               // product particle identifier lists
vector<string>   sigma_name [N_SPECIES];                     // identifier string for each process
vector<int>      target_to_species;                          // mapping target index to species index

// particle coordinates, containers aligned to multiple of 64 byte memory addresses to help vectorizations

const int        INIT_N_P = 100000;                          // initial number of particles for storage allocation (electrons / ions)
alignas(64) vector<double> x  [N_SPECIES];                   // coordinates of electrons (one spatial, three velocity components)
alignas(64) vector<double> vx [N_SPECIES];                   // coordinates of electrons (one spatial, three velocity components)
alignas(64) vector<double> vy [N_SPECIES];                   // coordinates of electrons (one spatial, three velocity components)
alignas(64) vector<double> vz [N_SPECIES];                   // coordinates of electrons (one spatial, three velocity components)

vector<bool>     coll_select;                                // helper vector for collision partner selection
vector<int>      coll_indexes [N_SPECIES];                   // storing indexes of colliding particles
vector<int>      out_indexes [N_SPECIES];                    // storing indexes of particles leaving

double           bias{};                                     // DC bias   
typedef array<double,N_G> xvector;                           // array for quantities defined at gird points
xvector          efield, pot;                                // electric field and potential
xvector          density [N_SPECIES];                        // particle densities
xvector          cumul_density [N_SPECIES];                  // cumulative densities
xvector          pfield,temp;                                // thermal power input and temperature field


xvector          target_density [N_TARGET];                  // target (background) species densities
double           max_target_density [N_TARGET];              // peak target (background) species densities

xvector          meta_source, meta_source_last;              // source distribution for metastables

typedef unsigned long long int Ullong;                       // compact name for 64 bit unsigned integer
Ullong   N_arrived  [N_SPECIES][N_SIDES] = {};               // counter for particles absorbed at the electrodes
Ullong   N_ele_emit [N_SPECIES][N_SIDES] = {};               // counter for emitted electrons at the electrodes
Ullong   N_ele_refl [N_SIDES]            = {};               // counter for reflected electrons at the electrodes

// electron energy probability function

array<double,N_EEPF> eepf;                                   // time integrated EEPF in the center of the plasma

// ion flux-energy distributions
array<int,N_FED> fed_pow [N_SPECIES];                        // FED at the powered electrode
array<int,N_FED> fed_gnd [N_SPECIES];                        // FED at the grounded electrode
double           mean_energy_pow [N_SPECIES] = {0.0};        // mean energy at the powered electrode
double           mean_energy_gnd [N_SPECIES] = {0.0};        // mean ion energy at the grounded electrode

// spatio-temporal (XT) distributions

const int N_BIN  = 20;                                       // number of time steps binned for the XT distributions
const int N_XT   = N_T / N_BIN;                              // number of spatial bins for the XT distributions
using xt_distr = array<double,N_G*N_XT>;                     // array for XT distributions (decimal numbers)
xt_distr  pot_xt;                                            // XT distribution of the potential
xt_distr  efield_xt;                                         // XT distribution of the electric field
xt_distr  process_xt[SAVE_XT_TOT];                           // XT distributions for selected collision rates
xt_distr  n_xt[N_SPECIES];                                   // XT distributions for densities
xt_distr  u_xt[N_SPECIES];                                   // XT distributions for mean velocities
xt_distr  j_xt[N_SPECIES];                                   // XT distributions for current
xt_distr  power_xt[N_SPECIES];                               // XT distributions for power absorption
xt_distr  meane_xt[N_SPECIES];                               // XT distributions for mean energies
xt_distr  pfield_xt;                                         // XT distributions for thermal power input field

// For Boltzmann term analysis

xt_distr  T_par_xt[N_SPECIES];                               // XT distribution for parallel SPECIES temperature --> actually pressure(!!!) up to save_xt where divided by n_xt[SPECIES]
xt_distr  Pi_C_xt[N_SPECIES];                                // XT distribution for SPECIES momentum loss   
xt_distr  term,term1,term2;

const double MIN_X = 0.45 * L;
const double MAX_X = 0.55 * L;

// event counters

array<vector<Ullong>,N_SPECIES> event_counter;               // arrays for counting individual collision events 

// Verboncoeour type solution of the Poisson equation

const double alpha_0 =  3.0 * RESISTANCE / (2.0 * DT_E) + INV_C;
const double alpha_1 = -2.0 * RESISTANCE / DT_E;
const double alpha_2 =  RESISTANCE / (2.0 * DT_E);
double       old_sigma;
double       old_charge;
double       old_old_charge;
double       qpow;
double       Kprime;

double   mean_energy_accu_center = 0.0;                      // mean electron energy accumulator in the center of the gap
Ullong   N_center_mean_energy    = 0;                        // mean electron energy counter in the center of the gap  
Ullong   N_coll [N_SPECIES]      = {0};                      // counter for collisions

double   Time;                                               // total simulated time (from the beginning of the simulation)
int      cycle,no_of_cycles,cycles_done;                     // current cycle and total cycles in the run, cycles completed (from the beginning of the simulation)

FILE     *datafile;                                          // used for saving data

bool     measurement_mode;                                   // flag that controls measurements and data saving

bool     DEBUG_MODE = false;
double   mean_temp  = T_WALL;                                // mean gas temperature
double   FA_E_THRESHOLD  = FA_E_THRESHOLD_FACTOR*mean_temp;  // fast atom thershold energy calculation
Ullong   FA_ADDED_ION        = 0;
Ullong   FA_ADDED_FA         = 0;
Ullong   FA_REMOVED          = 0;
Ullong   ENERGY_ADDED_ION    = 0;
Ullong   ENERGY_ADDED_FA     = 0;
Ullong   ENERGY_REMOVED_ION  = 0;
Ullong   ENERGY_REMOVED_FA   = 0;
Ullong   ENERGY_ADDED_FA_REM = 0;

//---------------------------------------------------------------------------//
// C++ Mersenne Twister 19937 generator                                      //
// R01(MTgen) will genarate uniform distribution over [0,1) interval         //
// RMB(MTgen) will generate Maxwell-Boltzmann distribution (of gas atoms)    //
//---------------------------------------------------------------------------//

std::random_device rd{}; 
std::mt19937 MTgen(rd());
std::uniform_real_distribution<> R01(0.0, 1.0);
std::normal_distribution<> RMB(0.0, 1.0);

double RNDveloc (const double mass, const double temp = mean_temp){
    return RMB(MTgen) * sqrt(K_BOLTZMANN * temp / mass);
}

//----------------------------------------------------------------------------//
//  Import cross sections from pre-processed data files                       //
//----------------------------------------------------------------------------//

void read_cross_sections(void){
    vector<float> buffer;
    char          st [CS_STRING_LENGTH];
    ifstream      file;
    double        cs_de, p_Eth;
    float         p_buff;
    int           cs_n, cs_sets, p_type, p_target, Np, p_product[4];
    string        csname, interpolation_type, comments, filename;
    bool          found;


    csname.reserve(CS_STRING_LENGTH);
    interpolation_type.reserve(CS_STRING_LENGTH);
    comments.reserve(CS_STRING_LENGTH);
    cout << ">> PICit! Loading cross-sections for " << N_SPECIES << " species" << endl;

    for (int species=0; species<N_SPECIES; species++){
        cout << "   Loading cross-sections from file " << CS_FILES[species] << endl;

        file.exceptions(ifstream::failbit | ifstream::badbit);
        try { file.open(string(CS_FILES[species]), ios::binary); }
        catch (system_error& e) { cerr << e.code().message() << endl; }

        // loading and printing common information
        file.read(st, CS_STRING_LENGTH); csname = st;      
        file.read(st, CS_STRING_LENGTH); interpolation_type = st;
        file.read(st, CS_STRING_LENGTH); comments = st;
        file.read((char*) &cs_de,   sizeof(double));
        file.read((char*) &cs_n,    sizeof(int));
        file.read((char*) &cs_sets, sizeof(int));
        if(cs_n  != CS_RANGES) cout << "Warning: Cross-section data lentgh inconsistent !!!" << endl;
        if(cs_de != DE_CS)     cout << "Warning: Cross-section data resolution inconsistent !!!" << endl;
        cout << "     Cross-section name:  " << csname << endl;
        cout << "     Interpolation type:  " << interpolation_type << endl;
        cout << "     Additional comments: " << comments << endl;
        cout << "     Energy resolution:   " << cs_de << " eV" << endl;
        cout << "     Total data sets:     " << cs_sets << endl;
        cout << "     Data points per set: " << cs_n << endl;
        cs_names[species] = comments;

        // allocating storage space
        N_CS[species] = cs_sets;                           
        sigma_type[species].reserve(N_CS[species]);
        sigma_target[species].reserve(N_CS[species]);
        sigma_threshold[species].reserve(N_CS[species]);
        sigma_name[species].reserve(N_CS[species]);
        sigma[species].clear();
        sigma[species].reserve(N_CS[species]*CS_RANGES);
        buffer.clear();
        buffer.reserve(N_CS[species]*CS_RANGES);

        // loading the individual CS sets
        for (int set=0; set<N_CS[species]; set++){         
            file.read((char*) &p_type, sizeof(int));
            file.read((char*) &p_target, sizeof(int));
            file.read((char*) &p_Eth,  sizeof(double));
            file.read((char*) &p_product[0], 4*sizeof(int));
            file.read(st, CS_STRING_LENGTH);
            cout << "    " << right << setw(3) << set << ". > " << st << " / type: " << p_type 
                 << " / target: " << MATERIALS[p_target].name << " / E_th: " << scientific <<  p_Eth;
            int i = 0; while ((p_product[i] >= 0) && (i<4)){ cout << " / product(" << i << "): " << MATERIALS[p_product[i]].name; i++; }
            cout << endl;
            sigma_type[species].push_back(p_type);
            sigma_target[species].push_back(p_target);
            sigma_threshold[species].push_back(p_Eth);
            sigma_name[species].push_back(string(st));
            sigma_product[species][0].push_back(p_product[0]);
            sigma_product[species][1].push_back(p_product[1]);
            sigma_product[species][2].push_back(p_product[2]);
            sigma_product[species][3].push_back(p_product[3]);
            Np = 0; for(i=0; i<4; i++){ if (p_product[i] >= 0) Np++; }
            sigma_Nproduct[species].push_back(Np);
            for (i=0; i<CS_RANGES; i++){
                file.read((char*) &p_buff, sizeof(float));
                buffer.push_back(p_buff);
            } 
        }
        file.close();

        // changing data order for better memory coalescence 
        for (int i=0; i<CS_RANGES; i++){                   
            for (int set=0; set<N_CS[species]; set++){
                sigma[species].push_back( buffer.at(set*CS_RANGES + i) );
            }
        }

        // mapping target and product species to local particle indexes 
        for (int set=0; set<N_CS[species]; set++){
            found = false;
            for (int t{0}; t<N_TARGET; t++) 
                if (sigma_target[species].at(set) == TARGET_KIND[t]){
                    sigma_target[species].at(set) = t;
                    found = true;
                }
            if (!found) cout << "PICit! ERROR: target species mismatch with file " << CS_FILES[species] <<", as "<< MATERIALS[sigma_target[species].at(set)].name <<" is not found among targets!"<< endl;
            for (int j{0}; j<sigma_Nproduct[species].at(set); j++){
                found = false;
                for (int n{0}; n<N_SPECIES; n++) 
                    if (sigma_product[species][j].at(set) == SPECIES_KIND[n]){
                        sigma_product[species][j].at(set) = n;
                        found = true;
                    }
                if (!found) cout << "PICit! warning: product species mismatch with file " << CS_FILES[species] <<", as "<<MATERIALS[ sigma_product[species][j].at(set)].name <<" is not found among species!"<< endl;
            }
        }

        // initializing event counters
        event_counter.at(species).clear();
        for (int set=0; set<N_CS[species]; set++) event_counter.at(species).push_back(0);

        cout << "     Selected for XT-analysis: ";
        if (SAVE_XT_NUM[species]){
            for (int i{0}; i<SAVE_XT_NUM[species]; i++) { cout << SAVE_XT_PROCESS[species][i] << ((i<(SAVE_XT_NUM[species]-1)) ? ",  " : ""); }
        } else { cout << "none"; }
        cout << endl;
    }
    // mappint target to species indexes
    target_to_species.resize(N_TARGET);
    for (int t{0}; t<N_TARGET; t++){
        found = false;
        for (int n{0}; n<N_SPECIES; n++) 
            if (TARGET_KIND[t] == SPECIES_KIND[n]){
                target_to_species.at(t) = n;
                found = true;
                cout << "   Mappint target type "<<MATERIALS[ TARGET_KIND[t] ].name <<" (No. "<< t <<") to species kind "<< MATERIALS[ SPECIES_KIND[n] ].name <<" (No. "<< n <<")" << endl;
            }
        if (!found) target_to_species.at(t) = -1;
    }
}

//---------------------------------------------------------------------//
// find upper limit of collision frequencies for unit density          //
//---------------------------------------------------------------------//

double max_coll_freq(int species, int target){
    double e,v,nu,nu_max;
    nu_max = 0;
    for(int i{0}; i<CS_RANGES; i++){
        e  = i * DE_CS;
        //v  = sqrt(2.0 * e / FACTOR_MQ(species));
        v  = sqrt(2.0 * e * EV_TO_J / MASS(species));
        nu = v * sigma_tot[species][target].at(i);
        if (nu > nu_max) {nu_max = nu;}
    }
    return nu_max;
}

void calc_max_coll_probability(void){
    // calculate max_coll_probability for null-collisions
    double tmp [N_SPECIES] = {0.0};
    for (int t{0}; t<N_TARGET; t++){    
        max_target_density[t] = *max_element(target_density[t].begin(), target_density[t].end());
        for (int n{0}; n<N_SPECIES; n++){
            tmp[n] += max_coll_freq0[n][t] * DT[n] * max_target_density[t];
        }
    }
    for (int n{0}; n<N_SPECIES; n++){
        max_coll_probability[n] = 1.0 - exp(-tmp[n]);
    }
 
}

//----------------------------------------------------------------------//
//  calculation of total cross sections for all species                 //
//----------------------------------------------------------------------//

void calc_total_cross_sections(void){
    float cs_sum;
    for(int n{0}; n<N_SPECIES; n++){
        for(int t{0}; t<N_TARGET; t++){
            sigma_tot[n][t].reserve(CS_RANGES);
            for(int i=0; i<CS_RANGES; i++){
                cs_sum = 0.0;
                for (int j{0}; j<N_CS[n]; j++){
                    if (sigma_target[n].at(j) == t) cs_sum += sigma[n].at(i*N_CS[n]+j); 
                }
                sigma_tot[n][t].push_back(cs_sum);   // total microscopic cross section
            }
            max_coll_freq0[n][t] = max_coll_freq(n,t);
        }
    }

    cout << ">> PICit! interaction matrix (No. of cross-processes):" << endl;
    const char Separator    = ' ';
    const int  Width        = 13;
    int        processes[N_SPECIES][N_TARGET];

    for (int n{0}; n<N_SPECIES; n++)
        for (int t{0}; t<N_TARGET; t++){
            processes[n][t] = 0;
            for (int j{0}; j<N_CS[n]; j++) if (sigma_target[n].at(j) == t) processes[n][t]++;
        }

    cout << "   " << right << setw(Width) << setfill(Separator) << "proj \\ targ" << " | ";
    for (int t{0}; t<N_TARGET; t++) cout << right << setw(Width) << setfill(Separator) << MATERIALS[TARGET_KIND[t]].name << " | ";
    cout << endl;
    cout << "   "  << right << setw(Width) << setfill('-') << "-" << "---";
    for (int t{0}; t<N_TARGET; t++) cout << right << setw(Width+3) << setfill('-') << "-";
    cout << endl;
    for (int n{0}; n<N_SPECIES; n++){
        cout << "   "  << right << setw(Width) << setfill(Separator) << MATERIALS[SPECIES_KIND[n]].name << " | ";
        for (int t{0}; t<N_TARGET; t++) cout << right << setw(Width) << setfill(Separator) << processes[n][t] << " | ";
        cout << endl;
    }
    cout << "   "  << right << setw(Width) << setfill('-') << "-" << "---";
    for (int t{0}; t<N_TARGET; t++) cout << right << setw(Width+3) << setfill('-') << "-";
    cout << endl;

}

//----------------------------------------------------------------------//
//  test of cross sections for all species                              //
//----------------------------------------------------------------------//

void test_cross_sections(const int n ){
    string filename="cross_sections_"+to_string(n)+".dat";
    ofstream file(filename);       // cross sections saved in data file: cross_sections_NUMBEROFSPECIES.dat

    for(int i=0; i<CS_RANGES; i++){
        file << scientific << setw(14) << i*DE_CS;
        for(int j=0; j<N_TARGET; j++)  file << " " << scientific << setw(14) << sigma_tot[n][j].at(i);
        for(int j=0; j<N_CS[n];  j++)  file << " " << scientific << setw(14) << sigma[n].at(i*N_CS[n]+j);
        file << endl;
    }
    file.close();
    string python_com="python csplot.py "+filename+" "+to_string(N_CS[n])+" " + to_string(N_TARGET);
    auto a=system(python_com.c_str()); (void)a;
}

//---------------------------------------------------------------------//
// initialize vectors for data storage & XTs if meauserent is set      //
//---------------------------------------------------------------------//

void init_vectors(void){

    fill(eepf.begin(), eepf.end(), 0.0);
    fill(meta_source.begin(), meta_source.end(), 0.0);
    coll_select.resize(2*INIT_N_P);

    // initialize also XT-s if measurement mode is on
    if(measurement_mode){
        fill(pot_xt.begin(), pot_xt.end(), 0.0);
        fill(efield_xt.begin(), efield_xt.end(), 0.0);
        fill(pfield_xt.begin(), pfield_xt.end(), 0.0);
        for (int i{0}; i<SAVE_XT_TOT; i++)
            fill(process_xt[i].begin(), process_xt[i].end(), 0.0);
    }

    for(size_t n=0; n<N_SPECIES; n++){
        x[n].reserve(INIT_N_P);
        vx[n].reserve(INIT_N_P);
        vy[n].reserve(INIT_N_P);
        vz[n].reserve(INIT_N_P);
        out_indexes[n].reserve(INIT_N_P/10);
        coll_indexes[n].reserve(INIT_N_P/10);

        if(measurement_mode){
            fill(n_xt[n].begin(), n_xt[n].end(), 0.0);
            fill(j_xt[n].begin(), j_xt[n].end(), 0.0);
            fill(u_xt[n].begin(), u_xt[n].end(), 0.0);
            fill(power_xt[n].begin(), power_xt[n].end(), 0.0);
            fill(meane_xt[n].begin(), meane_xt[n].end(), 0.0);
            fill(T_par_xt[n].begin(), T_par_xt[n].end(), 0.0);
            fill(Pi_C_xt[n].begin(), Pi_C_xt[n].end(), 0.0);
        }

        fill(cumul_density[n].begin(), cumul_density[n].end(), 0.0);
        fill(fed_pow[n].begin(), fed_pow[n].end(), 0);
        fill(fed_gnd[n].begin(), fed_gnd[n].end(), 0);
    }
}

//----------------------------------------------------------------------//
// initialization of the simulation by placing a given number of        //
// electrons and ions at random positions between the electrodes        //
//----------------------------------------------------------------------//

void init_seed(int nseed, bool forseeding = false){
    cout<<">> PICit! Seeding particles..."<<endl;
    for(int n=0; n<N_SPECIES; n++){
        if(n==AR_FAST) break;   //break from loop if fast atom TODO
        for (int i=0; i<nseed; i++){
            x[n].push_back( L * R01(MTgen) );                                  // initial random position
            vx[n].push_back(0.0); vy[n].push_back(0.0); vz[n].push_back(0.0);  // initial velocity components
        }
        if(forseeding) break;   //break from loop after seeding electrons (n=0)
    }
}


//----------------------------------------------------------------------//
// Specific collision processes                                         //
//----------------------------------------------------------------------//

void isotropic_scattering(const int species_index, const double &E_threshold, double &eta, double &chi, double &g ){
    if(E_threshold != 0.0){
        double energy = 0.5 * MASS(species_index) * g * g;      // projectile energy
        energy = fabs(energy - E_threshold * EV_TO_J);          // subtract energy loss for excitation
        g = sqrt(2.0 * energy / MASS(species_index)); 
    }  
    eta = TWO_PI * R01(MTgen);
    chi = acos(1.0 - 2.0 * R01(MTgen));
}

void ionization_Opal(const int species_index, const int product_index, const double E_threshold, const double Opal_w, double &eta, double &chi, double &g, const double xe,
 const double ct, const double st, const double cp, const double sp, const double F2, const double wx, const double wy, const double wz){

    double energy = 0.5 * MASS(species_index) * g * g;      // projectile energy
    energy = fabs(energy - E_threshold * EV_TO_J);          // subtract energy loss of ionization
    double e_ej = Opal_w * tan(R01(MTgen) * atan(0.5*energy*J_TO_EV / Opal_w)) * EV_TO_J; // energy of the ejected particle
    double e_sc = fabs(energy - e_ej);                      // energy of scattered projectile after the collision
    g    = sqrt(2.0 * e_sc / MASS(species_index));          // relative velocity of scattered projectile
    double g2 = sqrt(2.0 * e_ej / MASS(species_index));     // relative velocity of ejected particle
    chi  = acos(sqrt(e_sc / energy));                       // scattering angle for scattered projectile
    double chi2 = acos(sqrt(e_ej / energy));                // scattering angle for ejected particles
    eta  = TWO_PI * R01(MTgen);                             // azimuthal angle for scattered projectile
    double eta2 = eta + PI;                                 // azimuthal angle for ejected particle
    double sc   = sin(chi2);
    double cc   = cos(chi2);
    double se   = sin(eta2);
    double ce   = cos(eta2);
    double gx   = g2 * (ct * cc - st * sc * ce);
    double gy   = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    double gz   = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    x[ELE].push_back( xe );                                 // add new electron
    vx[ELE].push_back( wx + F2 * gx );
    vy[ELE].push_back( wy + F2 * gy );
    vz[ELE].push_back( wz + F2 * gz );

    double ratio = WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[product_index]; // determine how many ions need to be produced

    double rmodr{}, rmodl{};
    size_t p{};
    rmodr  = xe * INV_DX;
    p      = static_cast<size_t>(rmodr);
    rmodr  -= floor(rmodr);                // right-side remainder
    rmodl  = 1.0-rmodr;                    // left-side remainder
    double r = R01(MTgen);
    double temperature = temp.at(p)*rmodl+temp.at(p+1)*rmodr;
    while(r < ratio){
        x[product_index].push_back( xe );                                                    // add new ion
        vx[product_index].push_back( RNDveloc(MASS(product_index), temperature) );           // velocity is sampled from background thermal distribution
        vy[product_index].push_back( RNDveloc(MASS(product_index), temperature) ); 
        vz[product_index].push_back( RNDveloc(MASS(product_index), temperature) ); 
        --ratio;
    }
}

void ionization_equalshare(const int species_index, const int product_index, const double E_threshold, double &eta, double &chi, double &g, const double xe,
 const double ct, const double st, const double cp, const double sp, const double F2, const double wx, const double wy, const double wz){

    double energy = 0.5 * MASS(species_index) * g * g;      // projectile energy
    energy = fabs(energy - E_threshold * EV_TO_J);
    double e_ej = 0.5*energy;
    double e_sc = e_ej;
    g    = sqrt(2.0 * e_sc / MASS(species_index));          // relative velocity of scattered projectile
    chi  = acos(1.0-2.0*R01(MTgen));                        // scattering angle for scattered projectile   
    double chi2 = acos(1.0-2.0*R01(MTgen));                 // scattering angle for ejected particle
    eta  = TWO_PI * R01(MTgen);                             // azimuthal angle for scattered projectile
    double eta2 = TWO_PI * R01(MTgen);                      // azimuthal angle for ejected particle
    double sc   = sin(chi2);
    double cc   = cos(chi2);
    double se   = sin(eta2);
    double ce   = cos(eta2);
    double gx   = g * (ct * cc - st * sc * ce);
    double gy   = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    double gz   = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    x[ELE].push_back( xe );                                 // add new electron
    vx[ELE].push_back( wx + F2 * gx );
    vy[ELE].push_back( wy + F2 * gy );
    vz[ELE].push_back( wz + F2 * gz );

    double ratio = WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[product_index]; // determine how many ions need to be produced
    double r = R01(MTgen);                                  // determine how many ions need to be produced
    
    double rmodr{}, rmodl{};
    size_t p{};
    rmodr  = xe * INV_DX;
    p      = static_cast<size_t>(rmodr);
    rmodr  -= floor(rmodr);                // right-side remainder
    rmodl  = 1.0-rmodr;                    // left-side remainder
    double temperature = temp.at(p)*rmodl+temp.at(p+1)*rmodr;
    while(r < ratio){
        x[product_index].push_back( xe );                    // add new ion
        vx[product_index].push_back( RNDveloc(MASS(product_index), temperature) );           // velocity is sampled from background thermal distribution
        vy[product_index].push_back( RNDveloc(MASS(product_index), temperature) ); 
        vz[product_index].push_back( RNDveloc(MASS(product_index), temperature) ); 
        --ratio;
    }
}

// @@@
void isotropic_scattering_fa(const int species_index, const double E_threshold, double &eta, double &chi, double &g, const double xe,
 const double ct, const double st, const double cp, const double sp, const double F1, const double wx, const double wy, const double wz, const double energy_1_t){

    if(E_threshold != 0.0){
        double energy = 0.5 * MASS(species_index) * g * g;      // projectile energy
        energy = fabs(energy - E_threshold * EV_TO_J);          // subtract energy loss for excitation
        g = sqrt(2.0 * energy / MASS(species_index)); 
    }  
    eta = TWO_PI * R01(MTgen);
    chi = acos(1.0 - 2.0 * R01(MTgen));

    double sc   = sin(chi);
    double cc   = cos(chi);
    double se   = sin(eta);
    double ce   = cos(eta);
    double gx   = g * (ct * cc - st * sc * ce);
    double gy   = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    double gz   = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    double vx_2_t = wx - F1 * gx;
    double vy_2_t = wy - F1 * gy;
    double vz_2_t = wz - F1 * gz;

    double energy_2_t = 0.5 *  MASS(AR_FAST) * (vx_2_t * vx_2_t + vy_2_t * vy_2_t + vz_2_t * vz_2_t);


    double rmodr{}, rmodl{};
    size_t p{};
    //cout << "Target " << energy_2_t*J_TO_EV << endl;

    if(energy_2_t > FA_E_THRESHOLD){ // add new fast atom and remove its energy from P field
        bool toGenerate = (species_index == AR_FAST) ? true : R01(MTgen) < WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[AR_FAST];
        if (toGenerate){ // TODO AR SPECIFIC
            x[AR_FAST].push_back( xe );
            vx[AR_FAST].push_back( vx_2_t ); 
            vy[AR_FAST].push_back( vy_2_t );
            vz[AR_FAST].push_back( vz_2_t );

            rmodr  = xe * INV_DX;
            p      = static_cast<size_t>(rmodr);
            rmodr  -= floor(rmodr);                  // right-side remainder
            rmodl  = 1.0 - rmodr;                    // left-side remainder
            double dp = energy_1_t/(ELECTRODE_AREA*DX)*WEIGHT(AR_FAST);
            //cout << "- " << dp << endl;
            pfield.at(p)   -= rmodl * dp;
            pfield.at(p+1) -= rmodr * dp;

            //if (species_index == AR_FAST) {FA_ADDED_FA++; ENERGY_REMOVED_FA++;}
            //else if (species_index == AR_ION) {FA_ADDED_ION++; ENERGY_REMOVED_ION++;}
        }     
    }
    else{ // atom not fast enough, add its extra energy to P field

        rmodr  = xe * INV_DX;
        p      = static_cast<size_t>(rmodr);
        rmodr  -= floor(rmodr);                // right-side remainder
        rmodl  = 1.0 - rmodr;                  // left-side remainder
        double dp = (energy_2_t-energy_1_t)/(ELECTRODE_AREA*DX)*WEIGHT(species_index);
        //cout << "+ " << dp << endl;
        pfield.at(p)   += rmodl * dp;
        pfield.at(p+1) += rmodr * dp;
        
        //if (species_index == AR_FAST) ENERGY_ADDED_FA++;
        //else if (species_index == AR_ION) ENERGY_ADDED_ION++;
    }
}

// @@@
void backward_scattering_fa(const int species_index, const double E_threshold, double &eta, double &chi, double &g, const double xe,
 const double ct, const double st, const double cp, const double sp, const double F1, const double wx, const double wy, const double wz, const double energy_1_t){

    chi  = PI;                                              // scattering angle for scattered projectile
    eta  = TWO_PI * R01(MTgen);                             // azimuthal angle for scattered projectile

    double sc   = sin(chi);
    double cc   = cos(chi);
    double se   = sin(eta);
    double ce   = cos(eta);
    double gx   = g * (ct * cc - st * sc * ce);
    double gy   = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    double gz   = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    double vx_2_t = wx - F1 * gx;
    double vy_2_t = wy - F1 * gy;
    double vz_2_t = wz - F1 * gz;

    double energy_2_t = 0.5 *  MASS(AR_FAST) * (vx_2_t * vx_2_t + vy_2_t * vy_2_t + vz_2_t * vz_2_t);     // TODO

    double rmodr{}, rmodl{};
    size_t p{};
    //cout << "Target " << energy_2_t*J_TO_EV << endl;
    //cout << "back weight_proj" << weight_proj << endl;

    if(energy_2_t > FA_E_THRESHOLD){ // add new fast atom and remove its energy from P field
        bool toGenerate = (species_index == AR_FAST) ? true : R01(MTgen) < WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[AR_FAST];
        if (toGenerate){ // TODO AR SPECIFIC
            x[AR_FAST].push_back( xe );
            vx[AR_FAST].push_back( vx_2_t ); 
            vy[AR_FAST].push_back( vy_2_t );
            vz[AR_FAST].push_back( vz_2_t );

            rmodr  = xe * INV_DX;
            p      = static_cast<size_t>(rmodr);
            rmodr  -= floor(rmodr);                  // right-side remainder
            rmodl  = 1.0 - rmodr;                    // left-side remainder
            double dp = energy_1_t/(ELECTRODE_AREA*DX)*WEIGHT(AR_FAST);
            //cout << "- " << dp << endl;
            pfield.at(p)   -= rmodl * dp;
            pfield.at(p+1) -= rmodr * dp;

            //if (species_index == AR_FAST) {FA_ADDED_FA++; ENERGY_REMOVED_FA++;}
            //else if (species_index == AR_ION) {FA_ADDED_ION++; ENERGY_REMOVED_ION++;}
        }     
    }
    else{ // atom not fast enough, add its extra energy to P field

        rmodr  = xe * INV_DX;
        p      = static_cast<size_t>(rmodr);
        rmodr  -= floor(rmodr);                // right-side remainder
        rmodl  = 1.0 - rmodr;                  // left-side remainder
        double dp = (energy_2_t-energy_1_t)/(ELECTRODE_AREA*DX)*WEIGHT(species_index);
        //cout << "+ " << dp << endl;
        pfield.at(p)   += rmodl * dp;
        pfield.at(p+1) += rmodr * dp;
        
        //if (species_index == AR_FAST) ENERGY_ADDED_FA++;
        //else if (species_index == AR_ION) ENERGY_ADDED_ION++;
    }
}


void backward_scattering(double &eta, double &chi){
      eta = TWO_PI * R01(MTgen);
      chi = PI;
}

void detachment(const int species_index, const int target_index, const int product_index, const int part_index, bool& is_lost,
const double &E_threshold, double &eta, double &chi, double &g,
const double &ct, const double &st, const double &cp, const double &sp, const double &F2, const double &wx, const double &wy, const double &wz){

    double prod_init_energy = 1.5*K_BOLTZMANN*mean_temp; // initial energy of product particle???
    double g1 = sqrt(2.0*prod_init_energy/MASS(product_index));
    
    double ratio = WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[product_index];  // determine how many targets need to vanish
    double r = R01(MTgen);     
    while(r<ratio){
        //create new product (electron)

        double chi2 = acos(1.0-2.0*R01(MTgen));                 // scattering angle for ejected particle
        double eta2 = TWO_PI * R01(MTgen);                      // azimuthal angle for ejected particle
        double sc   = sin(chi2);
        double cc   = cos(chi2);
        double se   = sin(eta2);
        double ce   = cos(eta2);
        double gx   = g1 * (ct * cc - st * sc * ce);
        double gy   = g1 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        double gz   = g1 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
        x[product_index].push_back( x[species_index].at(part_index) );                                 // add new electron
        vx[product_index].push_back( wx + F2 * gx );
        vy[product_index].push_back( wy + F2 * gy );
        vz[product_index].push_back( wz + F2 * gz );
        --ratio;
    }
    
    if(species_index==ELE){
        //Collide incoming electron
        double energy = 0.5 * MASS(species_index) * g * g;      // projectile energy
        energy = fabs(energy - E_threshold * EV_TO_J);
        g    = sqrt(2.0 * energy / MASS(species_index));          // relative velocity of scattered projectile
        chi  = acos(1.0-2.0*R01(MTgen));                        // scattering angle for scattered projectile   
        eta  = TWO_PI * R01(MTgen);                             // azimuthal angle for scattered projectile

        //delete target
        ratio = WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[target_index];  // determine how many targets need to vanish
        double r = R01(MTgen);  
        size_t p{}, cell_index=x[species_index].at(part_index)*INV_DX;
        while(r < ratio){
            //find first target in cell of projectile particle
            p = find_if(x[target_index].begin(),x[target_index].end(),[cell_index](const auto &x){
                return (x>=cell_index*DX && x<(cell_index+1)*DX);
            }) - x[target_index].begin();

            // delete target
            if(p!=x[target_index].size()){
                x[target_index].at(p)  =  x[target_index].back();  x[target_index].pop_back();
                vx[target_index].at(p) = vx[target_index].back(); vx[target_index].pop_back();
                vy[target_index].at(p) = vy[target_index].back(); vy[target_index].pop_back();
                vz[target_index].at(p) = vz[target_index].back(); vz[target_index].pop_back();
            }
            --ratio;
        }  

    }
    else{
        //delete species
        x[species_index].at(part_index)  =  x[species_index].back();  x[species_index].pop_back();
        vx[species_index].at(part_index) = vx[species_index].back(); vx[species_index].pop_back();
        vy[species_index].at(part_index) = vy[species_index].back(); vy[species_index].pop_back();
        vz[species_index].at(part_index) = vz[species_index].back(); vz[species_index].pop_back();
    }

    is_lost=true;
}


//Collisions in which species + target both vanish

void reactive_scattering(const int species_index, const int target_index, const int part_index, bool& is_lost){

    double ratio = WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[target_index];  // determine how many targets need to vanish
    double r = R01(MTgen);                  
    size_t p{0}, cell_index=x[species_index].at(part_index)*INV_DX;

    while(r < ratio){
        //find first target in cell of projectile particle
        p = find_if(x[target_index].begin(),x[target_index].end(),[cell_index](const auto &x){
            return (x>=cell_index*DX && x<(cell_index+1)*DX);
        }) - x[target_index].begin();

        // delete target
        if(p!=x[target_index].size()){
            x[target_index].at(p)  =  x[target_index].back();  x[target_index].pop_back();
            vx[target_index].at(p) = vx[target_index].back(); vx[target_index].pop_back();
            vy[target_index].at(p) = vy[target_index].back(); vy[target_index].pop_back();
            vz[target_index].at(p) = vz[target_index].back(); vz[target_index].pop_back();
        }
        --ratio;
    }

    // Delete species 
    x[species_index].at(part_index)  =  x[species_index].back();  x[species_index].pop_back();
    vx[species_index].at(part_index) = vx[species_index].back(); vx[species_index].pop_back();
    vy[species_index].at(part_index) = vy[species_index].back(); vy[species_index].pop_back();
    vz[species_index].at(part_index) = vz[species_index].back(); vz[species_index].pop_back();

    is_lost= true;
}


void diss_attachment(const int species_index, const int product_index, const int part_index, bool& is_lost){

    double xe= x[species_index].at(part_index);

    // Delete species 
    x[species_index].at(part_index)  =  x[species_index].back();  x[species_index].pop_back();
    vx[species_index].at(part_index) = vx[species_index].back(); vx[species_index].pop_back();
    vy[species_index].at(part_index) = vy[species_index].back(); vy[species_index].pop_back();
    vz[species_index].at(part_index) = vz[species_index].back(); vz[species_index].pop_back();

    //create product with given probability
    double ratio = WEIGHT_FACTORS[species_index]/WEIGHT_FACTORS[product_index];  // determine how many targets need to vanish
    double r = R01(MTgen);

    double rmodr{}, rmodl{};
    size_t p{};
    rmodr  = xe * INV_DX;
    p      = static_cast<size_t>(rmodr);
    rmodr  -= floor(rmodr);                // right-side remainder
    rmodl  = 1.0-rmodr;                    // left-side remainder
    double temperature = temp.at(p)*rmodl+temp.at(p+1)*rmodr;
    while(r<ratio){
        x[product_index].push_back( xe );                    // add new ion
        vx[product_index].push_back( RNDveloc(MASS(product_index), temperature) );           // velocity is sampled from background thermal distribution
        vy[product_index].push_back( RNDveloc(MASS(product_index), temperature) ); 
        vz[product_index].push_back( RNDveloc(MASS(product_index), temperature) ); 
        --ratio;
    }

    is_lost = true;
}

//----------------------------------------------------------------------//
// General collision routine                                            //
//----------------------------------------------------------------------//

void do_collision(const int species_index, const size_t part_index, int *e_index, 
                    double *vx_t, double *vy_t, double *vz_t, double *t_dens){


    // decide collision process
    int N_cs = N_CS[species_index];
    vector<double> t(N_cs,0);
    int target_index = sigma_target[species_index].at(0);
    t.front() = (double)sigma[species_index].at(e_index[target_index]*N_cs) * t_dens[target_index];
    for(int k{1}; k<N_cs; ++k){
        target_index = sigma_target[species_index].at(k);
        t.at(k) = t.at(k-1) + (double)sigma[species_index].at(e_index[target_index]*N_cs+k) * t_dens[target_index];
    }
    int coll_index = lower_bound(t.begin(),t.end(),float(R01(MTgen))*t.back()) - t.begin();
    target_index = sigma_target[species_index].at(coll_index);
    event_counter.at(species_index).at(coll_index)++;


    const double F1 = MASS(species_index)  / (MASS(species_index) + TARGET_MASS(target_index));
    const double F2 = 1.0 - F1;

    bool is_lost = false;   //for detachment 

    double xe   = x[species_index].at(part_index);
    double vx_p = vx[species_index].at(part_index);
    double vy_p = vy[species_index].at(part_index);
    double vz_p = vz[species_index].at(part_index);

    // COM frame for projectile (relative velocity)
    double gx = vx_p-vx_t[target_index];
    double gy = vy_p-vy_t[target_index];
    double gz = vz_p-vz_t[target_index];
    double g  = hypot(gy,gx,gz);

    // COM frame for projectile (COM velocity)
    double wx = F1*vx_p + F2*vx_t[target_index];
    double wy = F1*vy_p + F2*vy_t[target_index];
    double wz = F1*vz_p + F2*vz_t[target_index];

    // find Euler angles
    double theta{}, phi{};

    if (gx == 0) { theta = 0.5 * PI; } else { theta = atan2(sqrt(gy * gy + gz * gz),gx); }
    if (gy == 0) { if (gz > 0){ phi = 0.5 * PI; } else { phi = - 0.5 * PI; }
    } else { phi = atan2(gz, gy); }

    double st = sin(theta);
    double ct = cos(theta);
    double sp = sin(phi);
    double cp = cos(phi);

    // Angles to be decided during collision
    double eta{}, chi{};

    // call specific scattering model
    switch(sigma_type[species_index].at(coll_index)){
        case COLL_ISO: {
            isotropic_scattering(species_index, sigma_threshold[species_index].at(coll_index), eta, chi, g);
        } 
            break;
        case COLL_BACK: {backward_scattering(eta, chi);}
            break;
        case COLL_ANISO:
            break;
        case COLL_ION:
            if (species_index == ELE) { 
                int heavy_product = sigma_product[species_index][sigma_Nproduct[species_index].at(coll_index)-1].at(coll_index);
                // if (target_index == 1) cout << "e+Ne ION: " << heavy_product << endl;
                ionization_equalshare(species_index, heavy_product, sigma_threshold[species_index].at(coll_index), eta, chi, g, xe, ct, st, cp, sp, F2, wx, wy, wz);
            }
            //else { heavy ionization }
            break;
        case COLL_DISS: isotropic_scattering(species_index, sigma_threshold[species_index].at(coll_index), eta, chi, g);
            break;
        case COLL_DET:{
            int product = sigma_product[species_index][0].at(coll_index);
            detachment(species_index,target_to_species[target_index],product,part_index,is_lost,sigma_threshold[species_index].at(coll_index),eta,chi,g,ct,st,cp,sp,F2,wx,wy,wz);}
            break;
        case COLL_NEUT:  reactive_scattering(species_index,target_to_species[target_index],part_index,is_lost);
            break;
        case COLL_DISS_ATTACH: {
            int product = sigma_product[species_index][0].at(coll_index);
            diss_attachment(species_index,product,part_index,is_lost);}
            break;
        case COLL_DISS_EXC: isotropic_scattering(species_index, sigma_threshold[species_index].at(coll_index), eta, chi, g);
            break;
        case COLL_DISS_REC: reactive_scattering(species_index,target_to_species[target_index],part_index,is_lost);
            break;
        case COLL_ASS_DET: {
            int product = sigma_product[species_index][0].at(coll_index);
            diss_attachment(species_index,product,part_index,is_lost);}
            break;
        case COLL_ISO_FASTATOM: {
                double energy_t = 0.5 * MASS(AR_FAST) * (vx_t[target_index] * vx_t[target_index] + vy_t[target_index] * vy_t[target_index] + vz_t[target_index] * vz_t[target_index]);
                isotropic_scattering_fa(species_index, sigma_threshold[species_index].at(coll_index), eta, chi, g, xe, ct, st, cp, sp, F1, wx, wy, wz, energy_t);
            }
            break;
        case COLL_BACK_FASTATOM: {
                double energy_t = 0.5 * MASS(AR_FAST) * (vx_t[target_index] * vx_t[target_index] + vy_t[target_index] * vy_t[target_index] + vz_t[target_index] * vz_t[target_index]);
                backward_scattering_fa(species_index, sigma_threshold[species_index].at(coll_index), eta, chi, g, xe, ct, st, cp, sp, F1, wx, wy, wz, energy_t);
            }
            break;
        default: cout << "Warning: collision type could not be identified !!!" << endl;
    }

    if(!is_lost){
        double sc = sin(chi);
        double cc = cos(chi);
        double se = sin(eta);
        double ce = cos(eta);

        // compute new relative velocity:
        
        gx = g * (ct * cc - st * sc * ce);
        gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    
        // post-collision velocity of the colliding electron
        
        double vx_p_2 = wx + F2 * gx;
        double vy_p_2 = wy + F2 * gy;
        double vz_p_2 = wz + F2 * gz;

        vx[species_index].at(part_index) = vx_p_2;
        vy[species_index].at(part_index) = vy_p_2;
        vz[species_index].at(part_index) = vz_p_2;

        if(species_index == AR_FAST){    // if fast atom slowed down below threshold, remove it 
            double energy_p_2 = 0.5 *  MASS(AR_FAST) * (vx_p_2 * vx_p_2 + vy_p_2 * vy_p_2 + vz_p_2 * vz_p_2);     // TODO

            
            //double energy_p_1 = 0.5 * MASS(AR_FAST) * (vx_p * vx_p + vy_p * vy_p + vz_p * vz_p);
            //cout << "Projectile " << energy_p_2 - energy_p_1 << endl;

            if ( energy_p_2 < FA_E_THRESHOLD){
                x[AR_FAST].at(part_index)  =  x[AR_FAST].back();  x[AR_FAST].pop_back();
                vx[AR_FAST].at(part_index) = vx[AR_FAST].back(); vx[AR_FAST].pop_back();
                vy[AR_FAST].at(part_index) = vy[AR_FAST].back(); vy[AR_FAST].pop_back();
                vz[AR_FAST].at(part_index) = vz[AR_FAST].back(); vz[AR_FAST].pop_back();

                // add its energy back to P field
                double rmodr{}, rmodl{};
                size_t p{};

                rmodr  = xe * INV_DX;
                p      = static_cast<size_t>(rmodr);
                rmodr  -= floor(rmodr);                // right-side remainder
                rmodl  = 1.0-rmodr;                    // left-side remainder
                double dp = (energy_p_2)/(ELECTRODE_AREA*DX)*WEIGHT(AR_FAST);
                pfield.at(p)   += rmodl * dp;
                pfield.at(p+1) += rmodr * dp;

                //ENERGY_ADDED_FA_REM++;
                //FA_REMOVED++;
            }
        }
    }

}

//-----------------------------------------------------------------//
// solve Poisson equation (tridiagonal matrix algorithm)           //
//-----------------------------------------------------------------//

void solve_Poisson (const xvector &rho1, double Voltage, double conv_charge=0.0){
    constexpr double A =  1.0;
    constexpr double B = -2.0;
    constexpr double C =  1.0;
    constexpr double S =  1.0 / (2.0 * DX);
    constexpr double ALPHA = -DX * DX / EPSILON0;
    xvector      g, w, f;
    int          i;
    
    // apply potential to the electrodes - boundary conditions
    pot.back()  = 0.0;                        // potential at the grounded electrode
    
    // solve Poisson equation
    
    for(i=1; i<=N_G-2; i++) f.at(i) = ALPHA * rho1.at(i);

    if(alpha_0!=0){
        const double B0 = -1.0-DX/(alpha_0*EPSILON0*ELECTRODE_AREA);
        Kprime = alpha_1 * old_charge + alpha_2 * old_old_charge;
        //   cout<<"qpow="<<qpow/WEIGHT/E_CHARGE<<endl;
        f.front() = ALPHA * (0.5 * rho1.front() + INV_DX * old_sigma + INV_DX / ELECTRODE_AREA * (conv_charge - old_charge + (Voltage - Kprime) / alpha_0));
        w.front() = C/B0;
        g.front() = f.front()/B0;
        for(i=1; i<=N_G-2; i++){
            w.at(i) = C / (B - A * w.at(i-1));
            g.at(i) = (f.at(i) - A * g.at(i-1)) / (B - A * w.at(i-1));
        }
        pot.at(N_G-2) = g.at(N_G-2);
        for (i=N_G-3; i>=0; i--) pot.at(i) = g.at(i) - w.at(i) * pot.at(i+1);  
    }
    else{
        pot.front()=Voltage;
        f.at(1) -= pot.at(0);
        f.at(N_G-2) -= pot.at(N_G-1);
        w.at(1) = C/B;
        g.at(1) = f.at(1)/B;
        for(i=2; i<=N_G-2; i++){
            w.at(i) = C / (B - A * w.at(i-1));
            g.at(i) = (f.at(i) - A * g.at(i-1)) / (B - A * w.at(i-1));
        }
        pot.at(N_G-2) = g.at(N_G-2);
        for (i=N_G-3; i>0; i--) pot.at(i) = g.at(i) - w.at(i) * pot.at(i+1);  
    } 

    // compute electric field
    for(i=1; i<=N_G-2; i++) efield.at(i) = (pot.at(i-1) - pot.at(i+1)) * S;      // electric field at the grid points between the electrodes
    efield.front() = (pot.at(0)     - pot.at(1))     * INV_DX - rho1.at(0)     * DX / (2.0 * EPSILON0);   // powered electrode
    efield.back()  = (pot.at(N_G-2) - pot.at(N_G-1)) * INV_DX + rho1.at(N_G-1) * DX / (2.0 * EPSILON0);   // grounded electrode
}



//-----------------------------------------------------------------//
// solve heat equation (Thomas algorithm)                      @@@ //
//-----------------------------------------------------------------//

void solve_heat (xvector pfield1, double tt, double Tg_0, double Tg_L, int iter){

    if(iter >= MAX_ITER){
        printf(">> eduPIC: Error: Maximum number of iteration reached while solving the heat equation.\n");
        exit(0);
    }

    printf("Temp in %f %f, iteration: %d\n", Tg_0, Tg_L, iter);
    const double A =  1.0;
    const double B = -2.0;
    const double C =  1.0;
    const double S = 1.0 / (2.0 * DX);
    const double ALPHA = -DX * DX / (KAPPA*tt);
    const double LAMBDADX = LAMBDA/DX;
    xvector temp2;
    xvector      g, w, f;
    int          i;

    for(i=0; i<N_G; i++){
        temp2.at(i) = temp.at(i);
    }
    
    temp2.front()     = Tg_0;
    temp2.back()      = Tg_L; 
    

    // solve Heat equation
    
    for(i=1; i<=N_G-2; i++) f.at(i) = ALPHA * pfield1.at(i);
    f.at(1) -= temp2[0];
    f.at(N_G-2) -= temp2[N_G-1];
    w.at(1) = C/B;
    g.at(1) = f.at(1)/B;
    for(i=2; i<=N_G-2; i++){
        w.at(i) = C / (B - A * w.at(i-1));
        g.at(i) = (f.at(i) - A * g.at(i-1)) / (B - A * w.at(i-1));
    }
    temp2.at(N_G-2) = g.at(N_G-2);
    for (i=N_G-3; i>0; i--) temp2.at(i) = g.at(i) - w.at(i) * temp2.at(i+1);            // potential at the grid points between the electrodes

    double Tg_0_calc = (temp2.at(1)-temp2.front())*LAMBDADX + T_WALL;
    double Tg_L_calc = (temp2.at(N_G-2)-temp2.back())*LAMBDADX + T_WALL;
    
    if(abs(Tg_0_calc-temp2.front())/Tg_0_calc > 1e-6 || abs(Tg_L_calc-temp2.back())/Tg_L_calc > 1e-6){
        solve_heat(pfield1, tt, Tg_0_calc, Tg_L_calc, iter+1);
    }
    else{
        for(i=0; i<N_G; i++){
            temp.at(i) = temp2.at(i);
        }
        printf("Final temp %f %f\n", temp.front(), temp.back());
    }
    
}

void mean_temp_calc(){
    mean_temp = 0;
    for (int k=0; k<N_G; k++){
        mean_temp += temp.at(k);
    }
    mean_temp /= N_G;
    FA_E_THRESHOLD  = 9.0*3.0/2.0*K_BOLTZMANN*mean_temp;       // energy threshold
}

//---------------------------------------------------------------------//
// calculate particle densitiy distributions from particle positions   //
//---------------------------------------------------------------------//

void calc_particle_density(int species){
    double   rmod;
    int      p;

    fill(density[species].begin(), density[species].end(), 0.0);          
    for(const auto& xx: x[species]){
        rmod  = xx * INV_DX;
        p     = static_cast<int>(rmod);
        rmod -= floor(rmod);
        density[species].at(p)   += (1.0-rmod) * FACTOR_D(species);             // assignment of the particles's to the grid points which enclose the particle
        density[species].at(p+1) += rmod * FACTOR_D(species);
    }
    density[species].front() *= 2.0;
    density[species].back()  *= 2.0;
    transform( cumul_density[species].begin(), cumul_density[species].end(), 
               density[species].begin(), cumul_density[species].begin(), plus<double>() );    
}


//---------------------------------------------------------------------//
// move particles in electric field                                    //
//---------------------------------------------------------------------//

void (*move_particles)(int);           // pointer to the proper move_particle function

void move_particles_no_B(int species){
    const int  Vsize = 8;
    int        Np    = x[species].size();
    int        VNp   = Np & (-Vsize);  // AND-ing with -Vsize will round down to the nearest lower multiple of Vsize. This works only if Vsize is a power of 2
    int        p, q;

    out_indexes[species].clear();

    // vectorized loop 
    Vec8d   Vx, VEl, VEr, Vvx, Vex, Vr, Vrr, Vout;
    Vec8i   Vp;
    double *Px  = x[species].data();
    double *Pvx = vx[species].data();
    int     pp[Vsize] __attribute__((aligned(32)));
    int     Nout;
    for(int k{0}; k<VNp; k+=Vsize){ 
        Vx.load(Px);    
        Vvx.load(Pvx);    
        VEr = Vx * INV_DX;
        VEl = truncate(VEr);
        Vr  = VEr - VEl;
        Vrr = 1.0 - Vr;
        Vp  = truncate_to_int32(VEl);
        Vp.store(pp);        
        for (int i{0}; i<Vsize; i++){
            VEl.insert(i, efield[pp[i]]);
            VEr.insert(i, efield[pp[i]+1]);           
        }
        Vex  = Vrr * VEl;
        Vex  = mul_add(Vr, VEr, Vex);
        Vvx += Vex * FACTOR_P(species);
        Vx  += Vvx * DT[species];
        Vvx.store(Pvx);
        Vx.store(Px);
        Px  += Vsize;
        Pvx += Vsize;
        // do boundary check
        Vout = select(Vx <= 0.0,  static_cast<Vec8d>(1.0), 0.0);
        Vout = if_add(Vx >= L,   Vout, 1.0);
        Vp   = truncate_to_int32(Vout);
        Nout = horizontal_count(Vp == 1);
        while (Nout) {
            q = horizontal_find_first(Vp == 1);
            Nout--;
            Vp.insert(q, 0);
            out_indexes[species].push_back(k+q);
        }
    }

    // scalar loop for remaining particles
    double e_x;
    double rmod;
    for(int k{VNp}; k<Np; k++){                             
        rmod  = x[species].at(k) * INV_DX;
        p     = static_cast<int>(rmod);
        rmod -= floor(rmod);
        e_x  = (1.0-rmod) * efield.at(p) + rmod * efield.at(p+1);
        vx[species].at(k) += e_x * FACTOR_P(species);
        x[species].at(k)  += vx[species].at(k) * DT[species];
        if ((x[species].at(k)<=0.0) || (x[species].at(k)>=L)) out_indexes[species].push_back(k);
    }
}

void move_particles_perpendicular_B(int species){
    const double OMEGA_C = B_FIELD / FACTOR_MQ(species);
    const double  CC = cos(OMEGA_C * DT[species]);
    const double  SS = sin(OMEGA_C * DT[species]);
    const int  Vsize = 8;
    int        Np    = x[species].size();
    int        VNp   = Np & (-Vsize);  // AND-ing with -Vsize will round down to the nearest lower multiple of Vsize. This works only if Vsize is a power of 2
    int        p, q;

    out_indexes[species].clear();

    // vectorized loop 
    Vec8d   Vx, VEl, VEr, Vvx, Vvxb, Vvy, Vex, Vr, Vrr, Vout;
    Vec8i   Vp;
    double *Px  = x[species].data();
    double *Pvx = vx[species].data();
    double *Pvy = vy[species].data();
    int     pp[Vsize] __attribute__((aligned(32)));
    int     Nout;
    for(int k{0}; k<VNp; k+=Vsize){ 
        Vx.load(Px);    
        Vvx.load(Pvx);    
        Vvy.load(Pvy);    
        VEr = Vx * INV_DX;
        VEl = truncate(VEr);
        Vr  = VEr - VEl;
        Vrr = 1.0 - Vr;
        Vp  = truncate_to_int32(VEl);
        Vp.store(pp);        
        for (int i{0}; i<Vsize; i++){
            VEl.insert(i, efield[pp[i]]);
            VEr.insert(i, efield[pp[i]+1]);           
        }
        Vex  = Vrr * VEl;
        Vex  = mul_add(Vr, VEr, Vex);
        Vvxb = mul_add(0.5*FACTOR_P(species), Vex, Vvx);
        Vvx  = CC * Vvxb;
        Vvx  = mul_add(SS, Vvy, Vvx);
        Vvx  = mul_add(0.5*FACTOR_P(species), Vex, Vvx);
        Vvy  = CC * Vvy;
        Vvy  = mul_add(-1.0*SS, Vvxb, Vvy);
        Vx   = mul_add(DT[species], Vvx, Vx);
        Vvx.store(Pvx);
        Vvy.store(Pvy);
        Vx.store(Px);
        Pvx += Vsize;
        Pvy += Vsize;
        Px  += Vsize;
        // do boundary check
        Vout = select(Vx <= 0.0,  static_cast<Vec8d>(1.0), 0.0);
        Vout = if_add(Vx >= L,   Vout, 1.0);
        Vp   = truncate_to_int32(Vout);
        Nout = horizontal_count(Vp == 1);
        while (Nout) {
            q = horizontal_find_first(Vp == 1);
            Nout--;
            Vp.insert(q, 0);
            out_indexes[species].push_back(k+q);
        }
    }

    // scalar loop for remaining particles
    double e_x, vx_b;
    double rmod;
    for(int k{VNp}; k<Np; k++){                             
        rmod  = x[species].at(k) * INV_DX;
        p     = static_cast<int>(rmod);
        rmod -= floor(rmod);
        e_x  = (1.0-rmod) * efield.at(p) + rmod * efield.at(p+1);
        vx_b = vx[species].at(k) + 0.5 * FACTOR_P(species) * e_x;
        vx[species].at(k)  =  vx_b * CC + vy[species].at(k) * SS + 0.5 * FACTOR_P(species) * e_x;
        vy[species].at(k)  = -vx_b * SS + vy[species].at(k) * CC;
        x[species].at(k)  += vx[species].at(k) * DT[species];
        if ((x[species].at(k)<=0.0) || (x[species].at(k)>=L)) out_indexes[species].push_back(k);
    }
}


//---------------------------------------------------------------------//
// boundary checking and handling function variants                    //
//---------------------------------------------------------------------//

void check_boundary_constant(int species){     // variant using constant R_ele, S_ele and gamma coefficients
    bool   out;
    double pos, R, v_sqr, energy; 
    size_t energy_index, side;
    int    p;
    
    for(int k=(out_indexes[species].size()-1); k>=0; k--){
        
        p   = out_indexes[species].at(k);
        
        pos = x[species].at(p);

        out = true;
        if (pos <= 0.0) {    // the paritcle is out at the powered electrode
            side = 0; 
            N_arrived[species][side]++;
            if(alpha_0 != 0.0) qpow += CHARGE(species)/E_CHARGE;               // charge arriving the powered electrode
            if (measurement_mode) {
                v_sqr  = SQR(vx[species].at(p)) + SQR(vy[species].at(p)) + SQR(vz[species].at(p));
                energy = 0.5 * v_sqr * FACTOR_MQ(species);
                energy_index = static_cast<size_t>(energy / DE_FED);
                if (energy_index < N_FED) { fed_pow[species].at(energy_index)++; }  // save FED at the grounded electrode
            }
            if (species == AR_FAST) x[AR_FAST].at(p) *= -1;
        }     
        if (pos >= L  ) {    // the paritcle is out at the grounded electrode
            side = 1; 
            N_arrived[species][side]++;
            if (measurement_mode) {
                v_sqr  = SQR(vx[species].at(p)) + SQR(vy[species].at(p)) + SQR(vz[species].at(p));
                energy = 0.5 * v_sqr * FACTOR_MQ(species);
                energy_index = static_cast<size_t>(energy / DE_FED);
                if (energy_index < N_FED) { fed_gnd[species].at(energy_index)++; }  // save FED at the grounded electrode
            }
            if (species == AR_FAST) x[AR_FAST].at(p) = 2*L-x[AR_FAST].at(p);
        }
        if (species==ELE) {    // manage electron induced processes
            R = R01(MTgen);
            if (R < SURF_R_ELE[side]) {                                       // elastic electron reflection
                if ((side==0) && (alpha_0!=0)) qpow += 1.0;                   // for Qconv reflected electron should not count.
                x[ELE].at(p) = side*L + SURF_NORMAL[side]*DX*0.001;
                vx[ELE].at(p) = SURF_NORMAL[side] * fabs(vx[ELE].at(p));
                N_ele_refl[side]++;
                out = false;
            } else if (R < SURF_E_EMISSION[ELE][side]+SURF_R_ELE[side]) {      // secondary electron emission
                if ((side==0) && (alpha_0!=0)) qpow += 1.0;                    // for Qconv emitted electron should count oppositely.
                x[ELE].push_back( side*L + SURF_NORMAL[side]*DX*0.001 );       // add new electron
                vx[ELE].push_back( SURF_NORMAL[side] * fabs(RNDveloc(E_MASS)) );
                vy[ELE].push_back( RNDveloc(E_MASS) );
                vz[ELE].push_back( RNDveloc(E_MASS) );
                N_ele_emit[ELE][side]++;
            }
        }
        else if (species == AR_FAST) {
            v_sqr  = SQR(vx[species].at(p)) + SQR(vy[species].at(p)) + SQR(vz[species].at(p));
            double E_i = 0.5 * MASS(AR_FAST) * v_sqr;
            double E_r = (1-ALPHA)*E_i+E_WALL*ALPHA;
            double c0 = sqrt(E_r/E_i);
            //cout << "c0 " << c0 << endl;
            vx[AR_FAST].at(p) *= -c0;
            vy[AR_FAST].at(p) *= c0;
            vz[AR_FAST].at(p) *= c0;
            out = false;
        }
        else {    // manage ion induced processes
            if (SURF_E_EMISSION[species][side] > 0.0){
                double P_emit = SURF_E_EMISSION[species][side] * WEIGHT_FACTORS[species]/WEIGHT_FACTORS[ELE];
                R = R01(MTgen);
                while(R <= P_emit){
                    if ((side==0) && (alpha_0!=0)) qpow += 1.0;                // for Qconv emitted electron should count oppositely.
                    x[ELE].push_back( side*L + SURF_NORMAL[side]*DX*0.001 );   // add new electron
                    vx[ELE].push_back( SURF_NORMAL[side] * fabs(RNDveloc(E_MASS)) );
                    vy[ELE].push_back( RNDveloc(E_MASS) );
                    vz[ELE].push_back( RNDveloc(E_MASS) );
                    N_ele_emit[species][side]++;
                    --P_emit;
                }
            }
        }
        if (out) {    // if still out: remove the particle
            x[species].at(p)  =  x[species].back();  x[species].pop_back();
            vx[species].at(p) = vx[species].back(); vx[species].pop_back();
            vy[species].at(p) = vy[species].back(); vy[species].pop_back();
            vz[species].at(p) = vz[species].back(); vz[species].pop_back();
        }
    }
}


//---------------------------------------------------------------------//
// measurement routine                                                 //
//---------------------------------------------------------------------//

void do_measurement(const int species, int t){
    size_t t_index = t/N_BIN;
    // first do the common part -- only for species == ELE
    if(species==ELE){
        // potential - XT
        transform(next(pot_xt.begin(), t_index*N_G),next(pot_xt.begin(), (t_index+1)*N_G),pot.begin(),next(pot_xt.begin(), t_index*N_G),plus<double>()); 
        // electric field - XT
        transform(next(efield_xt.begin(), t_index*N_G),next(efield_xt.begin(), (t_index+1)*N_G),efield.begin(),next(efield_xt.begin(), t_index*N_G),plus<double>()); 
    }
    // density - XT
    transform(next(n_xt[species].begin(), t_index*N_G),next(n_xt[species].begin(), (t_index+1)*N_G),density[species].begin(),next(n_xt[species].begin(), t_index*N_G),plus<double>());
    
    double mean_v{}, rmodr{}, rmodl{}, e_x{}, velocity{}, v2{}, energy{}, rate{}, ratep{};
    size_t p{}, energy_index{};

    int    proc_index, XT_counter{0}; 
    for(int i{0}; i<species; i++) XT_counter += SAVE_XT_NUM[i];

    for(size_t k{0}; k<x[species].size(); ++k){
        // calculate X for XT, plus take the time shift from the leapfrog into account
        rmodr  = x[species].at(k) * INV_DX;
        p      = static_cast<size_t>(rmodr);
        rmodr  -= floor(rmodr);                // right-side remainder
        rmodl  = 1.0-rmodr;                    // left-side remainder
        e_x    = rmodl * efield.at(p) + rmodr * efield.at(p+1);
        mean_v = vx[species].at(k) - 0.5*e_x*FACTOR_P(species);
        // calculate current density -- mean velocity can be calculated from j(x,t)!
        j_xt[species].at(t_index*N_G+p)   += rmodl*mean_v;
        j_xt[species].at(t_index*N_G+p+1) += rmodr*mean_v;

        // calculate kinetic pressure, p_xx
        v2 = SQR(mean_v);
        T_par_xt[species].at(t_index*N_G+p)   += rmodl*v2;
        T_par_xt[species].at(t_index*N_G+p+1) += rmodr*v2;
        
        // calculate energy of particle
        v2     = SQR(mean_v) + SQR(vy[species].at(k)) + SQR(vz[species].at(k));
        energy = 0.5*MASS(species)*v2*J_TO_EV;
        

        // mean energy -- XT
        meane_xt[species].at(t_index*N_G+p)   += rmodl*energy;
        meane_xt[species].at(t_index*N_G+p+1) += rmodr*energy;
        
        // calculate momentum loss, Pi_C_xt 
        energy_index = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES-1);
        velocity     = sqrt(v2);
        rate         = 0.0;
        ratep        = 0.0;
        for (int t{0}; t<N_TARGET; t++){
            rate += static_cast<double>(sigma_tot[species][t].at(energy_index))*target_density[t].at(p);
            ratep += static_cast<double>(sigma_tot[species][t].at(energy_index))*target_density[t].at(p+1);
        } 
        rate        *= mean_v*velocity;
        ratep        *= mean_v*velocity;
        Pi_C_xt[species].at(t_index*N_G+p)   -= rmodl*rate;
        Pi_C_xt[species].at(t_index*N_G+p+1) -= rmodr*ratep;
        
        // process rate -- XT
        for (int i{0}; i<SAVE_XT_NUM[species]; i++){
            proc_index = SAVE_XT_PROCESS[species][i];
            if (energy > sigma_threshold[species].at(proc_index)){
                rate   = static_cast<double>(sigma[species].at(energy_index*N_CS[species]+proc_index));
                rate  *= velocity;
                process_xt[XT_counter+i].at(t_index*N_G+p)   += rmodl*rate*target_density[sigma_target[species].at(proc_index)].at(p);
                process_xt[XT_counter+i].at(t_index*N_G+p+1) += rmodr*rate*target_density[sigma_target[species].at(proc_index)].at(p+1);
            }
        }
        
        if(species==ELE){
            // EEPF and mean energy at center
            if ((MIN_X < x[species][k]) && (x[species][k] < MAX_X)){
                energy_index = static_cast<int>(energy / DE_EEPF);
                if (energy_index < N_EEPF) {eepf[energy_index] += 1.0;}
                mean_energy_accu_center += energy;
                N_center_mean_energy++;
            }
        }
    }
}

//---------------------------------------------------------------------//
// Accumulate source function for metastables (electron impact only)   //
//---------------------------------------------------------------------//

void acc_meta_source(void){

    const size_t Ne = x[ELE].size();
    double mean_v{}, rmodr{}, rmodl{}, e_x{}, v2{}, energy{}, rate{};
    size_t p{}, energy_index{};

    for(size_t k{0}; k<Ne; ++k){
        // calculate X plus take the time shift from the leapfrog into account
        rmodr  = x[ELE].at(k) * INV_DX;
        p      = static_cast<size_t>(rmodr);
        rmodr  -= floor(rmodr);                // right-side remainder
        rmodl  = 1.0-rmodr;                    // left-side remainder
        e_x    = rmodl * efield.at(p) + rmodr * efield.at(p+1);
        mean_v = vx[ELE].at(k) - 0.5*e_x*FACTOR_P(ELE);

        // calculate energy of particle
        v2           = SQR(mean_v) + SQR(vy[ELE].at(k)) + SQR(vz[ELE].at(k));
        energy       = 0.5*FACTOR_MQ(ELE)*v2;

        // process rate -- XT
        if (energy > sigma_threshold[ELE].at(META_SOURCE_PROCESS)){
            energy_index = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES-1);
            rate         = static_cast<double>(sigma[ELE].at(energy_index*N_CS[ELE]+META_SOURCE_PROCESS));
            rate        *= sqrt(v2);
            meta_source.at(p)   += rmodl*rate*target_density[sigma_target[ELE].at(META_SOURCE_PROCESS)].at(p);
            meta_source.at(p+1) += rmodr*rate*target_density[sigma_target[ELE].at(META_SOURCE_PROCESS)].at(p+1);
        }
    }
}

//---------------------------------------------------------------------//
// Voltage waveform functions to choose from                           //
//---------------------------------------------------------------------//

double Voltage_cosine(int time_step){
    double phase = TWO_PI * static_cast<double>(time_step) / static_cast<double>(N_T);
    //return VOLTAGE * cos(phase) + bias;
    return VOLTAGE * cos(phase);
}

double Voltage_peaks(int time_step, const size_t harmonics){
    double V{}, phase = TWO_PI * static_cast<double>(time_step) / static_cast<double>(N_T);
    double r = 1.0 / SQR(static_cast<double>(harmonics)+1);
    
    for(size_t k{0}; k<harmonics; ++k){
        V += (harmonics-k)*cos(phase*(k+1));
    }
    return 2*V*VOLTAGE*r+bias;
}

double Voltage_valleys(int time_step, const size_t harmonics){
    double V{}, phase= TWO_PI * static_cast<double>(time_step) / static_cast<double>(N_T);
    double r = 1.0 / SQR(static_cast<double>(harmonics)+1);

    for(size_t k{0}; k<harmonics; ++k){
        V += (k%2) ? (harmonics-k)*cos(phase*(k+1)+PI) :(harmonics-k)*cos(phase*(k+1));
    }
    return 2*V*VOLTAGE*r+bias;
}


//---------------------------------------------------------------------//
// DC bias determination                                               //
//---------------------------------------------------------------------//

void calculate_bias(double& bias){
    double bl{}, br{};

    double ediff_l = WEIGHT_FACTORS[ELE]*static_cast<double>(N_arrived[ELE][0]-N_ele_refl[0]);     // absorbed electrons at the powered electrode
    double ediff_r = WEIGHT_FACTORS[ELE]*static_cast<double>(N_arrived[ELE][1]-N_ele_refl[1]);     // absorbed electrons at the grounded electrode
    double ionsum_l{0.0}, ionsum_r{0.0};                // other species

    for(size_t k{0}; k<N_SPECIES; ++k){               // subtract emitted electrons
        ediff_l -= WEIGHT_FACTORS[ELE]*N_ele_emit[k][0];
        ediff_r -= WEIGHT_FACTORS[ELE]*N_ele_emit[k][1];
        if(k>0){
            ionsum_l += WEIGHT_FACTORS[k]*N_arrived[k][0]*CHARGE(k)/E_CHARGE;
            ionsum_r += WEIGHT_FACTORS[k]*N_arrived[k][1]*CHARGE(k)/E_CHARGE;
        }
    }
    
    if(ediff_l+ionsum_l > 0){ bl = (ionsum_l-ediff_l)/(ediff_l+ionsum_l); }
    else{ bl = 1.0; }
    if(ediff_r+ionsum_r > 0){ br = (ionsum_r-ediff_r)/(ediff_r+ionsum_r); }
    else{ br = 1.0; }

    bias += (bl-br)*0.5;
}


//---------------------------------------------------------------------------//
// TEMPORARY VARIABLES FOR O2!!!!!!!!!!!!!                                   //
//---------------------------------------------------------------------------//


// void calc_o2a_density(void){

//     const double radius = 0.25;   // experimental reactor radius
//     double lambda_sq,v_th,v_over_a,diff,simtime,plus,rec_k;
//     const double  surface_loss_coeff = 6.0e-3;  // O2a destruction, wall recombination coefcient

//     lambda_sq = 1.0 / (pow(PI/L,2.0) + pow(2.4/radius,2.0));
//     v_th = sqrt((8.0*K_BOLTZMANN*TEMPERATURE)/(PI*MATERIALS[O2_GAS].mass));
//     v_over_a = L / (2.0 * (1.0 + L/radius));
//     diff = 0.2e-4 * 1e5  / PRESSURE;
//     simtime = no_of_cycles * PERIOD;
//     plus = metacounter*WEIGHT(ELE)/ (ELECTRODE_AREA * L * simtime);
//     rec_k = lambda_sq/diff + v_over_a * 2.0 * (2.0-surface_loss_coeff) / v_th / surface_loss_coeff;  // reciprocal of rate
//     o2a_density = plus * rec_k;
//     o2a_density_ratio = o2a_density / GAS_DENSITY; 

//     FILE   * f;
//     f = fopen("o2a_calc.txt","w");
//     fprintf(f,"## O2a density calculation ##\n");
//     fprintf(f,"electrode radius   = %e\n",radius); 
//     fprintf(f,"O2a creation rate  = %e\n",plus); 
//     fprintf(f,"surface loss coeff = %e\n",surface_loss_coeff); 
//     fprintf(f,"lambda_sq / diff   = %e\n",lambda_sq/diff); 
//     fprintf(f,"1/loss             = %e\n",rec_k); 
//     fprintf(f,"O2a density        = %e\n",o2a_density); 
//     fprintf(f,"O2a density ratio  = %e\n",o2a_density_ratio);
//     fclose(f);
// }


//---------------------------------------------------------------------------//
// TEMPORARY VARIABLES FOR O2!!!!!!!!!!!!!                                   //
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------//
// simulation of one radiofrequency cycle                              //
//---------------------------------------------------------------------//

void do_one_cycle (void){
    
    int      p, rint, Np, energy_index[N_TARGET], Nc[N_SPECIES];
    double   g, g_sqr, vx_a[N_TARGET], vy_a[N_TARGET], vz_a[N_TARGET], energy, nu, p_coll, Ncr_mod, rmod, t_dens[N_TARGET];
    double   new_charge{0.0}, new_sigma{0.0}, gen_voltage{0.0};
    xvector  rho;

    for (int t{0}; t<N_T; t++){                                  // the RF period is divided into N_T equal time intervals (time step DT_E)
        Time += DT_E;                                           // update of the total simulated time
        
        // step 1: compute densities at grid points
        calc_particle_density(ELE);                             // electron density - computed in every time step
        if ((t % N_SUB) == 0) {                                 // ion density - computed in every N_SUB-th time steps (subcycling)
            for(int i{1}; i<N_SPECIES; ++i) calc_particle_density(i); 
        }

        // step 2: solve Poisson equation
        
        fill(rho.begin(),rho.end(),0.0);
        for(int n{0}; n<N_G; ++n){
            for(int k{0}; k<N_SPECIES; ++k) rho[n] += CHARGE(k) * density[k][n];
        }
        
        gen_voltage=Voltage_cosine(t);                          // generator voltage
        if(alpha_0!=0){
            qpow *= E_CHARGE * WEIGHT(ELE); 
            solve_Poisson(rho, gen_voltage, qpow);              // compute potential and electric field
            new_charge = (gen_voltage + pot.back() - pot.front() - Kprime) / alpha_0;   // new charge
            new_sigma  = old_sigma + (qpow + new_charge - old_charge) / ELECTRODE_AREA; // new surface charge at the powered electrode
            old_old_charge = old_charge;                                                // replace old by new
            old_charge = new_charge;                                                    // replace old by new
            old_sigma  = new_sigma;                                                     // replace old by new
            qpow = 0.0;                                                                 // clear accumulator for the next time step
        } else {
            solve_Poisson(rho, gen_voltage);                    // compute potential and electric field
        }
        
        // steps 3 & 4: move particles according to electric field interpolated to particle positions
        move_particles(ELE);                                    // move all electrons in every time steps   
        if ((t % N_SUB) == 0) {                                 // move all ions in every N_SUB-th time steps (subcycling)
            for(int k{1}; k<N_SPECIES; ++k) move_particles(k); 
        }

        // step 5: check boundaries
        check_boundary_constant(ELE);                           // check boundaries for all electrons in every time step
        if ((t % N_SUB) == 0) {                                 // check boundaries for all ions in every N_s-th time steps (subcycling)
            for(int k{1}; k<N_SPECIES; ++k) check_boundary_constant(k);
        }
        
        //Do measurements
        if(measurement_mode){
            for(int k{0}; k<N_SPECIES; ++k) do_measurement(k,t);

//            for (p=0; p<N_G; p++) {
//                pfield_xt[p][t_index] += pfield[p];
//            }
        }
        
        // accumulate metastable source
        if(META_T_INDEX){
           acc_meta_source();
        }

        // step 6: collisions
        // NULL-collision: select scattering particles
        for(int n{0}; n<N_SPECIES; n++){
            if ((n == ELE) || ((t % N_SUB) == 0)){
                Np       = x[n].size();
                Ncr_mod  = Np * max_coll_probability[n];
                Nc[n]    = static_cast<int>(Ncr_mod);
                //if (n == AR_FAST) cout << "max_coll_freq0[AR_FAST][0] " << max_coll_freq0[AR_FAST][0] << endl;
                Ncr_mod -= floor(Ncr_mod);
                if (R01(MTgen) < Ncr_mod) Nc[n]++;
                coll_indexes[n].resize(Nc[n]);
                if (coll_select.size() < Np) coll_select.resize(Np+1000);
                fill(coll_select.begin(), next(coll_select.begin(), Np), false);
                for(int i{0}; i<Nc[n]; i++){
                    do{ p = Np * R01(MTgen); } while (coll_select.at(p));
                    coll_select.at(p) = true;
                    coll_indexes[n].at(i) = p;
                }
            }
        }
        
        // electron collision
        for (int b{0}; b<N_TARGET; b++){ vx_a[b] = 0.0; vy_a[b] = 0.0; vz_a[b] = 0.0; } // in case of cold gas model
        for (int k{0}; k<Nc[ELE]; k++){                                   // loop over colliding electrons
            p = coll_indexes[ELE].at(k);
            if(p >= x[ELE].size()) continue;
            rmod  = x[ELE].at(p) * INV_DX;
            rint  = static_cast<int>(rmod);
            rmod -= floor(rmod);
            nu    = 0.0;
            double temperature = temp.at(rint)*(1.0-rmod)+temp.at(rint+1)*rmod;
            for (int b{0}; b<N_TARGET; b++){                              // loop over target species
                t_dens[b] = (1.0-rmod) * target_density[b].at(rint) + rmod * target_density[b].at(rint+1);
                if (WARM_GAS) {                                           // pick velocity components of a random target gas atom
                    vx_a[b] = RNDveloc(TARGET_MASS(b), temperature); vy_a[b] = RNDveloc(TARGET_MASS(b), temperature); vz_a[b] = RNDveloc(TARGET_MASS(b), temperature); 
                }  
                g_sqr  = SQR(vx[ELE].at(p)-vx_a[b]) + SQR(vy[ELE].at(p)-vy_a[b]) + SQR(vz[ELE].at(p)-vz_a[b]);
                g      = sqrt(g_sqr);
                energy = 0.5 * g_sqr * FACTOR_MQ(ELE);
                energy_index[b] = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES-1);
                nu    += sigma_tot[ELE][b].at(energy_index[b]) * g * t_dens[b];
            }
            p_coll = 1.0 - exp(-nu * DT[ELE]);                        // collision probability for electrons
            if (R01(MTgen) < p_coll/max_coll_probability[ELE]){       // REAL collision  
                do_collision(ELE,p,energy_index,vx_a,vy_a,vz_a,t_dens);
                N_coll[ELE]++;
            }
        }
        
        // ion collisions
        if ((t % N_SUB) == 0) {                                           // checking for occurrence of a collision for all ions in every N_SUB-th time steps (subcycling)
            for(int n{1}; n<N_SPECIES; ++n){                              // loop over projectile species
                for (int k{0}; k<Nc[n]; k++){                             // loop over colliding projectiles
                    p = coll_indexes[n].at(k);
                    if(p >= x[n].size()) continue;
                    rmod  = x[n].at(p) * INV_DX;
                    rint  = static_cast<int>(rmod);
                    rmod -= floor(rmod);
                    nu    = 0.0;
                    double temperature = temp.at(rint)*(1.0-rmod)+temp.at(rint+1)*rmod;
                    for (int b{0}; b<N_TARGET; b++){                          // loop over target species
                        t_dens[b] = (1.0-rmod) * target_density[b].at(rint) + rmod * target_density[b].at(rint+1);
                        vx_a[b]   = RNDveloc(TARGET_MASS(b), temperature);                 // pick velocity components of a random target particle
                        vy_a[b]   = RNDveloc(TARGET_MASS(b), temperature); 
                        vz_a[b]   = RNDveloc(TARGET_MASS(b), temperature);  
                        g_sqr     = SQR(vx[n].at(p)-vx_a[b]) + SQR(vy[n].at(p)-vy_a[b]) + SQR(vz[n].at(p)-vz_a[b]);
                        g         = sqrt(g_sqr);
                        energy    = 0.5 * g_sqr * MU(n,b) * J_TO_EV;
                        energy_index[b] = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES-1);
                        nu       += sigma_tot[n][b].at(energy_index[b]) * g * t_dens[b];
                    }
                    p_coll = 1.0 - exp(-nu * DT[n]);                          // collision probability for ions
                    //if (n == AR_FAST) cout << "AR FAST " << p_coll/max_coll_probability[n] << endl;
                    if (R01(MTgen) < p_coll/max_coll_probability[n]) {        // REAL collision takes place
                        do_collision(n,p,energy_index,vx_a,vy_a,vz_a,t_dens);
                        N_coll[n]++;
                    }
                }
            }
        }
        
        if ((t % 1000) == 0) {
            if (INV_C==0) {
                printf(" c = %8d  t = %8d  #e = %8zd  ", cycle, t, x[ELE].size());
                for (int k{1}; k<N_SPECIES; k++)  printf("#i[%d] = %8zd  ", k, x[k].size());
                printf(" V = %8.3f\n", pot[0]);
            } else {
                printf(" c = %8d  t = %8d  #e = %8zd  ", cycle, t, x[ELE].size());
                for (int k{1}; k<N_SPECIES; k++)  printf("#i[%d] = %8zd  ", k, x[k].size());
                printf(" V = %8.3f VC = %8.3f\n",pot[0], old_charge*INV_C);
            }
        }
    }
    
    fprintf(datafile,"%8d  ",cycle);
    for(int k{0}; k<N_SPECIES; ++k) fprintf(datafile, "%8zd  ",x[k].size());
    fprintf(datafile,"%E\n",bias); 
}

//---------------------------------------------------------------------//
// save particle coordinates                                           //
//---------------------------------------------------------------------//

void update_target_densities(const int act_cycle, const int total_cycle){

    // placeholder for more sophisticted algorithms
    for (int t{0}; t<N_TARGET; t++){     
        //put constant density according to target_raito
        //if(TARGET_RATIO[t]>0.0) fill(target_density[t].begin(), target_density[t].end(), GAS_DENSITY*TARGET_RATIO[t]);
        for (int k=0; k<N_G; k++){
                target_density[t].at(k) = PRESSURE / (K_BOLTZMANN * temp.at(k));
            }
        for(size_t k{1};k<N_SPECIES;++k){
            // if t is one of species, target_density = cumul_density for given species
            if(SPECIES_KIND[k]==TARGET_KIND[t]) transform(target_density[t].begin(),target_density[t].end(),cumul_density[k].begin(),target_density[t].begin(),
            [act_cycle](auto x, auto y){
                (void)x; auto c = 1.0 / static_cast<double>(act_cycle) / static_cast<double>(N_T);
                return y*N_SUB*c;
            });
        }
    } 

    // metastable density calculation [http://dx.doi.org/10.1088/0022-3727/49/18/185202]
    if(META_T_INDEX && total_cycle){
        if ((total_cycle % META_CYCLE_SKIP) == 0){
            constexpr double A =  1.0;
            constexpr double B = -2.0;
            constexpr double C =  1.0;
            const double Diff_meta   = META_DIFF_0 * pow(mean_temp/300.0,1.5) / PRESSURE;   // Diffusion coeffient
            const double v_avg       = sqrt(3.0*K_BOLTZMANN*mean_temp/TARGET_MASS(META_T_INDEX));     // avg speed
	        const double lambda_meta = 2.0 * Diff_meta / v_avg;                                         // mean free path
	        const double beta_meta   = lambda_meta * (1.0 + META_REFLECTION) / (1.0 - META_REFLECTION) / SQRT3;
            const double ALPHA       = -DX * DX / Diff_meta;
            const double n0_factor   = beta_meta / (DX + beta_meta);
            xvector   g, w, f;
            int       i, n=0;
            double    n_ratio;

            transform(meta_source.begin(), meta_source.end(), meta_source_last.begin(), [](auto r){return r * FACTOR_D(ELE) / static_cast<double>(META_CYCLE_SKIP * N_T); });
            fill(meta_source.begin(), meta_source.end(), 0.0);
            fill(target_density[META_T_INDEX].begin(), target_density[META_T_INDEX].end(), max(target_density[META_T_INDEX].front(),1e16)/n0_factor); 
            do{ n++;
                for(i=0; i<N_G; i++) f.at(i) = ALPHA * meta_source_last.at(i);
                target_density[META_T_INDEX].front() *= n0_factor / (target_density[META_T_INDEX].front() / target_density[META_T_INDEX].at(1));
                target_density[META_T_INDEX].back()  *= n0_factor / (target_density[META_T_INDEX].back() / target_density[META_T_INDEX].at(N_G-2));
                f.at(1) -= target_density[META_T_INDEX].front();
                f.at(N_G-2) -= target_density[META_T_INDEX].back();
                w.at(1) = C/B;
                g.at(1) = f.at(1)/B;
                for(i=2; i<=N_G-2; i++){
                    w.at(i) = C / (B - A * w.at(i-1));
                    g.at(i) = (f.at(i) - A * g.at(i-1)) / (B - A * w.at(i-1));
                }
                target_density[META_T_INDEX].at(N_G-2) = g.at(N_G-2);
                for (i=N_G-3; i>0; i--) target_density[META_T_INDEX].at(i) = g.at(i) - w.at(i) * target_density[META_T_INDEX].at(i+1);  

                n_ratio = target_density[META_T_INDEX].front() / target_density[META_T_INDEX].at(1);
                if ((n%50) == 0) printf("META (%3d): slope %e, n0 %e   (peak %e)\n", n, n_ratio, target_density[META_T_INDEX].front(), target_density[META_T_INDEX].at(N_G/2) );
            } while (fabs((n_ratio / n0_factor) - 1.0) > 1e-6);
        }
    }

    calc_max_coll_probability();
}

//---------------------------------------------------------------------//
// save particle coordinates                                           //
//---------------------------------------------------------------------//

void save_particle_data(void){
    double   d;
    size_t   Np;
    FILE   * f;
    
    f = fopen("picdata.bin","wb");
    fwrite(&Time,sizeof(double),1,f);
    fwrite(&old_charge,sizeof(double),1,f);
    fwrite(&old_old_charge,sizeof(double),1,f);
    fwrite(&old_sigma,sizeof(double),1,f);
    fwrite(&qpow,sizeof(double),1,f);
    d = (double)(cycles_done);
    fwrite(&d,sizeof(double),1,f);
    for (int n=0; n<N_SPECIES; n++){
        Np = x[n].size();
        d  = static_cast<double>(Np);
        fwrite(&d,sizeof(double),1,f);
        fwrite(x[n].data(), sizeof(double),Np,f);
        fwrite(vx[n].data(),sizeof(double),Np,f);
        fwrite(vy[n].data(),sizeof(double),Np,f);
        fwrite(vz[n].data(),sizeof(double),Np,f);
    }

    // save temperature vector
    fwrite(&N_G,sizeof(double),1,f);
    fwrite(temp.data(), sizeof(double),N_G,f);

    fclose(f);
    printf(">> PICit! data saved : %zd electrons %zd ions, %d cycles completed, time is %e [s]\n", x[ELE].size(), x[1].size(), cycles_done, Time);
    printf(">> PICit! data saved : old_Q = %e, Q = %e, SIGMA = %e, QCONV = %e\n",old_old_charge,old_charge,old_sigma,qpow);
}

//---------------------------------------------------------------------//
// load particle coordinates                                           //
//---------------------------------------------------------------------//

void load_particle_data(void){
    double   d;
    FILE   * f;
    size_t   Np;
    
    f = fopen("picdata.bin","rb");
    if (f==NULL) {printf(">> PICit! error: No particle data file found, try running initial cycle using argument '0'\n"); exit(0); }
    auto a=fread(&Time,sizeof(double),1,f);
    a = fread(&old_charge,sizeof(double),1,f);
    a = fread(&old_old_charge,sizeof(double),1,f);
    a = fread(&old_sigma,sizeof(double),1,f);
    a = fread(&qpow,sizeof(double),1,f);
    a = fread(&d,sizeof(double),1,f);
    cycles_done = static_cast<int>(d);
    for (int n=0; n<N_SPECIES; n++){
        a = fread(&d,sizeof(double),1,f);
        Np = static_cast<size_t>(d);
        x[n].resize(Np);
        a = fread(x[n].data(),  sizeof(double),Np,f);
        vx[n].resize(Np);
        a = fread(vx[n].data(), sizeof(double),Np,f);
        vy[n].resize(Np);
        a = fread(vy[n].data(), sizeof(double),Np,f);
        vz[n].resize(Np);
        a = fread(vz[n].data(), sizeof(double),Np,f);
    }

    // load temperature vector
    a = fread(&d,sizeof(double),1,f);
    a = fread(temp.data(),  sizeof(double),N_G,f);

    (void)a;
    fclose(f);
    printf(">> PICit: data loaded : %zd electrons %zd ions, %d cycles completed before, time is %e [s]\n", x[ELE].size(), x[1].size(), cycles_done, Time);
    printf(">> PICit! data loaded : old_Q = %e, Q = %e, SIGMA = %e, QCONV = %e\n",old_old_charge,old_charge,old_sigma,qpow);
}

//---------------------------------------------------------------------//
// save / load target densities and meta_source                        //
//---------------------------------------------------------------------//

void save_target_densities(void){
    int      Ng = N_G;
    FILE   * f;
    
    f = fopen("pictarget.bin","wb");
    fwrite(&Ng, sizeof(int), 1, f);
    for (int t=0; t<N_TARGET; t++) fwrite(target_density[t].data(), sizeof(double), Ng, f);
    fwrite(meta_source.data(), sizeof(double), Ng, f);
    fclose(f);
    printf(">> PICit! target densities saved \n");
}

void load_target_densities(void){
    int      Ng;
    FILE   * f;
    
    f = fopen("pictarget.bin","rb");
    if (f==NULL) {printf(">> PICit! error: No target density file found, try running initial cycle using argument '0'\n"); exit(0); }
    auto a=fread(&Ng, sizeof(int), 1, f);
    if (Ng == N_G){
        for (int t{0}; t<N_TARGET; t++) a = fread(target_density[t].data(),  sizeof(double), Ng, f);
        a = fread(meta_source.data(),  sizeof(double), Ng, f);
        printf(">> PICit! target densities loaded \n");
    } else {
        printf(">> PICit! error: target density grid size mismatch, try running initial cycle using argument '0') \n");  exit(0); 
        // TODO: implement interpolation between different grids
    }
    (void)a;
    fclose(f);

    calc_max_coll_probability();

    printf("      max densities: ");
    for (int t{0}; t<N_TARGET; t++) printf(" %10.2e (%s)  ", max_target_density[t], MATERIALS[TARGET_KIND[t]].name.data());
    printf("\n");
    printf("      max coll probabilities: \n");
    for (int n{0}; n<N_SPECIES; n++){
        printf("      %12s:", MATERIALS[SPECIES_KIND[n]].name.data());
        printf(" %12.2e ", max_coll_probability[n]);
        printf("\n");
    }
}


//---------------------------------------------------------------------//
// save density data                                                   //
//---------------------------------------------------------------------//

void save_density(void){
        FILE *f;
        
        f = fopen("density.dat","w");
        double c = 1.0 / static_cast<double>(no_of_cycles) / static_cast<double>(N_T);
        for(size_t m=0; m<N_G; m++){
            fprintf(f,"%8.5f  ", m*DX);
            fprintf(f,"%12e  ", cumul_density[ELE].at(m) * c);
            for(size_t k{1};k<N_SPECIES;++k){
                fprintf(f,"%12e  ", cumul_density[k].at(m) * N_SUB * c);
            }
            fprintf(f,"%12e  ", target_density[META_T_INDEX].at(m));
            fprintf(f,"%12e  ", meta_source_last.at(m));
            fprintf(f,"\n");
        }
        fclose(f);
    }


void save_xvector(xvector& vec, const string& fname){
        FILE    *f;
        size_t   i,j;
        
        f = fopen(fname.c_str(),"w");
        for (i=0; i<N_G; i++){
            fprintf(f,"%8.5f  ", i*DX);
            fprintf(f,"%e", vec.at(i));
            fprintf(f,"\n");
        }
        fclose(f);
        cout << "   Saving " << fname << "...\n";
    }


//---------------------------------------------------------------------//
// save EEPF data                                                      //
//---------------------------------------------------------------------//

void save_eepf(void) {
    FILE   *f;
    double h,energy;
    
    h = accumulate(eepf.begin(), eepf.end(), 0.0);
    h *= DE_EEPF;
    f = fopen("eepf.dat","w");
    for (size_t i=0; i<N_EEPF; i++) {
        energy = (i + 0.5) * DE_EEPF;
        fprintf(f,"%e  %e\n", energy, eepf.at(i) / h / sqrt(energy));
    }
    fclose(f);
}

//---------------------------------------------------------------------//
// save IFED data                                                      //
//---------------------------------------------------------------------//

void save_fed(size_t species) {
    FILE   *f;
    double h_pow, h_gnd, energy;
    
    h_pow = accumulate(fed_pow[species].begin(), fed_pow[species].end(), 0.0);
    h_gnd = accumulate(fed_gnd[species].begin(), fed_gnd[species].end(), 0.0);
    h_pow *= DE_FED;
    h_gnd *= DE_FED;
    mean_energy_pow[species] = 0.0;
    mean_energy_gnd[species] = 0.0;
    if (species==0){ f = fopen("efed.dat","w"); } else { f = fopen("ifed.dat","w"); }
    for (size_t i=0; i<N_FED; i++) {
        energy = (i + 0.5) * DE_FED;
        fprintf(f,"%6.2f %10.6f %10.6f\n", energy, static_cast<double>(fed_pow[species].at(i))/h_pow, static_cast<double>(fed_gnd[species].at(i))/h_gnd);
        mean_energy_pow[species] += energy * static_cast<double>(fed_pow[species].at(i)) / h_pow;
        mean_energy_gnd[species] += energy * static_cast<double>(fed_gnd[species].at(i)) / h_gnd;
    }
    fclose(f);
}

//--------------------------------------------------------------------//
// save XT data                                                       //
//--------------------------------------------------------------------//

void save_xt_1(xt_distr& distr, const string& fname, double norm=1.0) {
    FILE    *f;
    size_t   i,j;
    
    f = fopen(fname.c_str(),"w");
    for (i=0; i<N_G; i++){
        for (j=0; j<N_XT; j++){
            distr.at(j*N_G+i) *= norm;
            fprintf(f,"%e  ", distr.at(j*N_G+i));
        }
        fprintf(f,"\n");
    }
    fclose(f);
    cout << "   Saving " << fname << "...\n";
}


//---------------------------------------------------------------------//
// Numerical differentiation routine for the Boltzmann term analysis   //
//---------------------------------------------------------------------//

//coordinate: 0=t,1=x
void differentiate(const xt_distr& input, xt_distr& output, const int coordinate){
    double fact{};
    size_t xin,xfin,tin,tfin;

    fill(output.begin(),output.end(),0.0);

    if(coordinate==0){ fact=N_XT/(12.0*PERIOD); xin=0; xfin=N_G; tin=2; tfin=N_XT-2;
    // use 5-stencil differentiation
        for(size_t t{tin};t<tfin;++t){ 
            for(size_t x{xin};x<xfin;++x){
                output.at(t*N_G+x) = (input.at(x+(t-2)*N_G)-8.0*input.at(x+(t-1)*N_G)+8.0*input.at(x+(t+1)*N_G)-input.at(x+(t+2)*N_G))*fact;
            }
        }
    }
    else{ fact=INV_DX/12.0; xin=2; xfin=N_G-2; tin=0; tfin=N_XT;
        for(size_t t{tin};t<tfin;++t){ 
            for(size_t x{xin};x<xfin;++x){
                output.at(t*N_G+x) = (input.at(x-2+t*N_G)-8.0*input.at(x-1+t*N_G)+8.0*input.at(x+1+t*N_G)-input.at(x+2+t*N_G))*fact;
            }
        }
    }
}


//---------------------------------------------------------------------//
// Full Boltzmann term analysis for arbitrary species                  //
//---------------------------------------------------------------------//

void Boltzmann_term_analysis(const int species, string sp_indent, const double norm){
    string fname{};
    

    // first calculate T_par properly, i.e. normalize by n_xt and subtract u_xt[species]**2, T_par=<v^2>-u^2
    transform(T_par_xt[species].begin(),T_par_xt[species].end(),n_xt[species].begin(),T_par_xt[species].begin(),[species](auto x, auto y){y/=FACTOR_D(species); if(y>1.0){return x/y;}else{return x;}});
    transform(T_par_xt[species].begin(),T_par_xt[species].end(),u_xt[species].begin(),T_par_xt[species].begin(),[](auto x, auto y){return x-y*y;});
    

    // Start calculating terms: first electric field, save, then multiply by current density.

    // TERM1, i.e. inertia term, E_in=m/nq[d(nu)/dt+d(nu^2)/dx]
    // calculate d(nu)/dt in term1
    differentiate(j_xt[species],term1,0);
    // calculate nu^2 in term2
    transform(j_xt[species].begin(),j_xt[species].end(),u_xt[species].begin(),term2.begin(),[](auto x, auto y){return x*y;});
    // calculate d(nu^2)/dx in term
    differentiate(term2,term,1);
    // term1+term, then divide by n --> TERM1E in term
    transform(term.begin(),term.end(),term1.begin(),term.begin(),plus<double>());
    // transform(term.begin(),term.end(),n_xt[species].begin(),term.begin(),[](auto x, auto y){y/=FACTOR_D; if(y>0.0){if(y>=1.0){return x/y;}else{return x;}}else return 0.0;});
    transform(term.begin(),term.end(),n_xt[species].begin(),term.begin(),[species](auto x, auto y){y/=FACTOR_D(species); if(y>0.0){return x/y;} else return 0.0;});
    // calculate TERM1H in term1
    transform(term.begin(),term.end(),j_xt[species].begin(),term1.begin(),[](auto& x, auto& y){return x*y;});
    // save files
    fname="TERME1_"+sp_indent+"_xt.dat";
    save_xt_1(term,fname,MASS(species)/CHARGE(species));
    fname="TERMH1_"+sp_indent+"_xt.dat";
    save_xt_1(term1,fname,MASS(species)*norm*FACTOR_D(species));
    // Inertia term done----------------------------------------

    // TERM2, i.e. ambipolar field, E_am=(T/nq)*(dn/dx)
    // calculate dn/dx, then multiply by T/n
    differentiate(n_xt[species],term,1);
    transform(term.begin(),term.end(),T_par_xt[species].begin(),term.begin(),[](auto x, auto y){return x*y;});
    // NOTE: as n_xt~FACTOR_D, here we do not divide!
    transform(term.begin(),term.end(),n_xt[species].begin(),term.begin(),[](auto x, auto y){ if(y>0.0){return x/y;}else return 0.0;});
    // TERM2H to term1
    transform(term.begin(),term.end(),j_xt[species].begin(),term1.begin(),[](auto x, auto y){return x*y;});
    // saving files
    fname="TERME2_"+sp_indent+"_xt.dat";
    save_xt_1(term,fname,MASS(species)/CHARGE(species)); 
    fname="TERMH2_"+sp_indent+"_xt.dat";
    save_xt_1(term1,fname,MASS(species)*norm*FACTOR_D(species));
    // TERM2 done------------------------------------------------------

    // TERM3, i.e. gradient T term, E_gradT=1/q*dT/dx
    differentiate(T_par_xt[species],term,1);
    transform(term.begin(),term.end(),j_xt[species].begin(),term1.begin(),[](auto x, auto y){return x*y;});
    // saving files
    fname="TERME3_"+sp_indent+"_xt.dat";
    save_xt_1(term,fname,MASS(species)/CHARGE(species)); 
    fname="TERMH3_"+sp_indent+"_xt.dat";
    save_xt_1(term1,fname,MASS(species)*norm*FACTOR_D(species));
    // TERM3 done---------------------------------------------------------------------------------------------------

    // Calculating Ohmic terms, i.e. TERM4=E_Ohm=-Pi_c/nq
    transform(Pi_C_xt[species].begin(),Pi_C_xt[species].end(),n_xt[species].begin(),term.begin(),[species](auto x, auto y){y/=FACTOR_D(species); if(y>0.0){return x/y;}else return 0.0;});
    transform(term.begin(),term.end(),j_xt[species].begin(),term1.begin(),[](auto x, auto y){return x*y;});
    // saving files
    fname="TERME4_"+sp_indent+"_xt.dat";
    save_xt_1(term,fname,-MASS(species)/CHARGE(species));
    fname="TERMH4_"+sp_indent+"_xt.dat";
    save_xt_1(term1,fname,-MASS(species)*norm*FACTOR_D(species));
    // Ohmic term done-------------------------------------------------------------------------------------------------

    // save T_par by multiplying by m_e/e to get in units of eV!
    fname="T_par"+sp_indent+"_xt.dat";
    save_xt_1(T_par_xt[species],fname,MASS(species)/E_CHARGE);

    // save Pi_C -- multiply by m_e 
    fname="Pi_C"+sp_indent+"_xt.dat";
    save_xt_1(Pi_C_xt[species],fname,FACTOR_D(species)*MASS(species)*norm);
}


//---------------------------------------------------------------------//
// Sheath analysis                                                     //
//---------------------------------------------------------------------//

void Sheath_analysis(void){
   // const double dN_G  = static_cast<double>(N_G);
    const double dN_XT = static_cast<double>(N_XT);
    const int ION = 1;

    int    i,j,k;
    FILE  *f;

    cout << "   Saving sheath_edge.dat" << endl;
	f = fopen("sheath_edge.dat","wt");

    // find sheath edge based on Brinkmann criterium [Brinkmann R P 2007 J. Appl. Phys. 102 093303]
    array<double,N_XT>  sheath_edge[N_SIDES];
    array<double,N_G>   dens_negative;
    array<double,N_G>   dens_positive;
    array<double,N_G>   dens_summ;
	int                 pos_max, si, sj;
	double              sum1, sum2;
    for(k = 0; k < N_G;  k++){ dens_positive[k]=0.0; dens_negative[k]=0.0; }
	for(i = 0; i < N_XT; i++){
        for(j = 0; j < N_SPECIES; j++){
            if(CHARGE(j) > 0.0) for(k = 0; k < N_G; k++) dens_positive[k] += n_xt[j].at(i*N_G+k);
            if(CHARGE(j) < 0.0) for(k = 0; k < N_G; k++) dens_negative[k] += n_xt[j].at(i*N_G+k);
        }
        for(k = 0; k < N_G; k++) dens_summ[k] = dens_positive[k] + dens_negative[k];
        pos_max = distance( dens_summ.begin(), max_element(dens_summ.begin(), dens_summ.end()) );
		si = 0;	
        sj = pos_max-1;
		sum1 = 0.5*dens_negative.at(si);
		sum2 = dens_positive.at(sj) - dens_negative.at(sj);
		while(si < sj){
			si++;
			sum1 += dens_negative.at(si);
			while ((sum1 > sum2) && (sj>1)){
				sj--;
				sum2 += dens_positive.at(sj) - dens_negative.at(sj);
			}
		}
        sheath_edge[0].at(i) = 0.5*max((si+sj),0) * DX;
		si = N_G-1;	
        sj = pos_max;
		sum1 = 0.5*dens_negative.at(si);
		sum2 = dens_positive.at(sj) - dens_negative.at(sj);
		while(si > sj){
			si--;
			sum1 += dens_negative.at(si);
			while ((sum1 > sum2) && (sj<N_G-2)){
				sj++;
				sum2 += dens_positive.at(sj) - dens_negative.at(sj);
			}
		}
        sheath_edge[1].at(i) = 0.5*min((si+sj),2*N_G-2) * DX;
		fprintf(f, "%10.6f   %.4e   %.4e   %.4e\n", (i+0.5)/dN_XT, sheath_edge[0].at(i), sheath_edge[1].at(i), pot_xt.at(i*N_G) );
	}
	fclose(f);

    // process applied voltage vaweform
    array<double,N_XT>  Phi;                                                             // Phi(t) = applied voltage vaweform
    for (size_t i=0; i<Phi.size(); i++) Phi.at(i) = pot_xt.at(i*N_G);
    double t_Phi_max = distance( Phi.begin(), max_element(Phi.begin(), Phi.end()) );     // time of maximum applied voltage 
    double t_Phi_min = distance( Phi.begin(), min_element(Phi.begin(), Phi.end()) );     // time of minimum applied voltage
    double Phi_avg   = accumulate(Phi.begin(), Phi.end(), 0.0) / Phi.size();             // measured DC level
    double Phi_max   = Phi.at(t_Phi_max) - Phi_avg;                                      // maximum applied voltage
    double Phi_min   = Phi.at(t_Phi_min) - Phi_avg;                                      // minimum applied voltage
 
    double s_max[N_SIDES];                                          // maximum sheat extensions (at max/min applied voltage)
    s_max[0] = sheath_edge[0].at(t_Phi_min);
    s_max[1] = L - sheath_edge[1].at(t_Phi_max);

    // calculate floating potentials, densitites and charges
    // full time-evolution is made available for future refinement based on true extrema instead of VW turning points
    array<double,N_XT>  sheath_voltage[N_SIDES];                    // voltage drop across the sheath 
    array<double,N_XT>  sheath_ion_densitiy[N_SIDES];               // ion density in the sheath
    array<double,N_XT>  sheath_charge[N_SIDES];                     // uncompensated charges in the sheath
    double              ni_mean[N_SIDES] = {0.0, 0.0};              // mean ion density in the sheath
    double              floating_pot[N_SIDES];                      // floating potentials at the powered and grounded electrodes
    double              Vs_max[N_SIDES];                            // maximum respective sheath voltages
    double              Q_max[N_SIDES];                             // maximum uncompensated charges in the sheath
    double              rmod, rint, Vs;
    size_t              p, s;

    // left sheath
    s = static_cast<size_t>(s_max[0]*INV_DX);
    for (size_t i=0; i<N_XT; i++){         
        rmod = modf(sheath_edge[0].at(i) * INV_DX, &rint);
        p    = static_cast<size_t>(rint);
        Vs = ((1.0-rmod)*pot_xt.at(i*N_G + p) + rmod*pot_xt.at(i*N_G + p+1));
        sheath_voltage[0].at(i) = Phi.at(i) - Vs;
        sheath_charge[0].at(i)  = 0.0;
        for(size_t j=0; j<=p; j++) for(size_t k=0; k<N_SPECIES; k++) sheath_charge[0].at(i) += n_xt[k].at(i*N_G+j) * CHARGE(k)/FACTOR_D(k); 
        sheath_ion_densitiy[0].at(i) = 0.0;
        for(size_t j=0; j<=s; j++) for(size_t k=1; k<N_SPECIES; k++) sheath_ion_densitiy[0].at(i) += n_xt[k].at(i*N_G+j); 
        sheath_ion_densitiy[0].at(i) /= static_cast<double>(s+1);
        ni_mean[0] += sheath_ion_densitiy[0].at(i);
    }
    floating_pot[0] = sheath_voltage[0].at(t_Phi_max);
    Vs_max[0]       = sheath_voltage[0].at(t_Phi_min);
    ni_mean[0]     /= dN_XT;
    Q_max[0]        = sheath_charge[0].at(t_Phi_min);

    // right sheath
    s = (N_G-1) - static_cast<size_t>(s_max[1]*INV_DX);
    for (size_t i=0; i<N_XT; i++){       
        rmod = 1.0 - modf(sheath_edge[1].at(i) * INV_DX, &rint);
        p    = static_cast<size_t>(rint);
        Vs = ((1.0-rmod)*pot_xt.at(i*N_G + p) + rmod*pot_xt.at(i*N_G + p-1));
        sheath_voltage[1].at(i) = Vs - pot_xt.at(i*N_G+(N_G-1));
        sheath_charge[1].at(i) = 0.0;
        for(size_t j=N_G-1; j>=p; j--) for(size_t k=0; k<N_SPECIES; k++) sheath_charge[1].at(i) += n_xt[k].at(i*N_G+j) * CHARGE(k)/FACTOR_D(k);  
        sheath_ion_densitiy[1].at(i) = 0.0;
        for(size_t j=N_G-1; j>=s; j--) for(size_t k=1; k<N_SPECIES; k++) sheath_ion_densitiy[1].at(i) += n_xt[k].at(i*N_G+j); 
        sheath_ion_densitiy[1].at(i) /= N_G-s;
        ni_mean[1] += sheath_ion_densitiy[1].at(i);
    }
    floating_pot[1] = sheath_voltage[1].at(t_Phi_min);
    Vs_max[1]       = sheath_voltage[1].at(t_Phi_max);
    ni_mean[1]     /= dN_XT;
    Q_max[1]        = sheath_charge[1].at(t_Phi_max);

    // maximum and minimum bulk potentials
    double bulk_pot_max    = Phi_max - floating_pot[0] - Vs_max[1];
    double bulk_pot_min    = Phi_min - floating_pot[1] - Vs_max[0];

    // sheath integrals
    double I_pow = 0.0;
    for (size_t i=0; i<=static_cast<size_t>(s_max[0]*INV_DX); i++){
        I_pow += n_xt[ION].at(t_Phi_min*N_G+i)/ni_mean[0] * (i+0.5)*DX/s_max[0] * DX/s_max[0];
    }
    I_pow *= 2.0;
    double I_gnd = 0.0;
    for (size_t i=N_G-1; i>=static_cast<size_t>((L-s_max[0])*INV_DX); i--){
        I_gnd += n_xt[ION].at(t_Phi_max*N_G+i)/ni_mean[1] * (N_G-0.5-i)*DX/s_max[1] * DX/s_max[1];
    }
    I_gnd *= 2.0;

    // terms of the symmetry parameter
    double ni_mean_ratio   = ni_mean[0] / ni_mean[1];
    double Q_max_ratio_sqr = SQR( Q_max[1] / Q_max[0] );
    double I_ratio         = I_gnd/I_pow;

    // symmetry parameter
    double epsilon = ni_mean_ratio * I_ratio * Q_max_ratio_sqr;

    // DC bias model results
    double eta_vw = -(Phi_max + epsilon*Phi_min)/(1.0+epsilon);
    double eta_f  = (floating_pot[0] + epsilon*floating_pot[1])/(1.0+epsilon);
    double eta_b  = (bulk_pot_max + epsilon*bulk_pot_min)/(1.0+epsilon);;
    double eta    = eta_vw + eta_f + eta_b;


    cout << "Phi_max/min:       " << fixed << setprecision(4) << Phi_max << "/" << fixed << setprecision(4) << Phi_min << endl;
    cout << "s_m/L (pow/gnd):   " << fixed << setprecision(4) << s_max[0]/L << "/" << fixed << setprecision(4) << s_max[1]/L << endl;
    cout << "Phi_f (pow/gnd):   " << fixed << setprecision(4) << floating_pot[0] << "/" << fixed << setprecision(4) << floating_pot[1] << endl;
    cout << "Phi_b (max/min):   " << fixed << setprecision(4) << bulk_pot_max << "/" << fixed << setprecision(4) << bulk_pot_min << endl;
    cout << "Vs_max (pow/gnd):  " << fixed << setprecision(4) << Vs_max[0] << "/" << fixed << setprecision(4) << Vs_max[1] << endl;
    cout << "Q_max (pow/gnd):   " << scientific << setprecision(4) << Q_max[0] << "/" << scientific << setprecision(4) << Q_max[1] << endl;

    cout << "-------------------------" << endl;
    cout << "nominal bias: " << fixed << setprecision(4) << bias << endl;
    cout << "<Phi(t)>:     " << fixed << setprecision(4) << Phi_avg << endl;
    cout << "eta:          " << fixed << setprecision(4) << eta << endl;
    cout << "eta_vw:       " << fixed << setprecision(4) << eta_vw << endl;
    cout << "eta_f:        " << fixed << setprecision(4) << eta_f << endl;
    cout << "eta_b:        " << fixed << setprecision(4) << eta_b << endl;
    cout << "-------------------------" << endl;
    cout << "epsilon:      " << fixed << setprecision(4) << epsilon << endl;
    cout << "ni_ratio:     " << fixed << setprecision(4) << ni_mean_ratio << endl;
    cout << "Q2_ratio:     " << fixed << setprecision(4) << Q_max_ratio_sqr << endl;
    cout << "Is_ratio:     " << fixed << setprecision(4) << I_ratio  << endl;
    cout << "-------------------------" << endl;

}



//---------------------------------------------------------------------//
// Current analysis at given position + Fourier transform              //
//---------------------------------------------------------------------//

void Current_analysis(const size_t pos){

    cout<<"   Current analysis started..."<<endl;

    const size_t f_limit=51;
    array<double, f_limit> volt_re, volt_im, curr_re, curr_im;
    fill(volt_re.begin(),volt_re.end(),0.0);
    fill(volt_im.begin(),volt_im.end(),0.0);
    fill(curr_re.begin(),curr_re.end(),0.0);
    fill(curr_im.begin(),curr_im.end(),0.0);

    // Print currents at position pos to Current_densities.dat
    FILE* f = fopen("Current_densities.dat","w");
    double ion_current{}, tot_current{}, ph{};
    fprintf(f,"#   t/T       voltage          electron           ion       displacement       total\n");

    for(int i=0;i<N_XT;++i){
        for(size_t k{1};k<N_SPECIES;++k) ion_current+=j_xt[k].at(i*N_G+pos);
        tot_current=j_xt[ELE].at(i*N_G+pos)+ion_current+term.at(i*N_G+pos);

        fprintf(f,"%8.4f  %14e  %14e  %14e  %14e  %14e\n",
        static_cast<double>(i)/static_cast<double>(N_XT),
                pot_xt.at(i*N_G), j_xt[ELE].at(i*N_G+pos), ion_current, term.at(i*N_G+pos), tot_current);

    //Fourier transform

        for(size_t k{0};k<f_limit;++k){
            ph = static_cast<double>(k*i) * 2.0 * PI / N_XT;
            volt_re.at(k) += pot_xt.at(i*N_G)*cos(ph)/N_XT;
            volt_im.at(k) += pot_xt.at(i*N_G)*sin(ph)/N_XT;
            curr_re.at(k) += tot_current*cos(ph)/N_XT;
            curr_im.at(k) += tot_current*sin(ph)/N_XT;
        }
        ion_current=0.0;
        tot_current=0.0;
    }
    fclose(f);

    // Save Fourier transform to file
    double va{}, ia{}, imp{}, phase{};
    f = fopen("Harmonics.dat","w");
    fprintf(f,"#k       V[k]           j[k]             Z[k]          angle[k]        Z_Re[k]         Z_Im[k]\n");
    for(size_t k{0};k<f_limit;k++){
        va = sqrt(volt_re.at(k)*volt_re.at(k)+volt_im.at(k)*volt_im.at(k));
        ia = sqrt(curr_re.at(k)*curr_re.at(k)+curr_im.at(k)*curr_im.at(k));
        imp = va / ia;
        phase = atan2(curr_im.at(k),curr_re.at(k)) - atan2(volt_im.at(k),volt_re.at(k));
        fprintf(f,"%zd  %14e  %14e  %14e  %14e  %14e  %14e  \n",
                k, va * 2.0, ia * 2.0, imp, phase, imp*cos(phase), imp*sin(phase));
    }
    fclose(f);

    cout<<"   Current analysis done..."<<endl;
}

//---------------------------------------------------------------------//
// Call analysis and save XT-s                                         //
//---------------------------------------------------------------------//

void save_all_xt(void){
    string fname;
    const double INV_COUNTS_PER_BIN = static_cast<double>(N_XT) / static_cast<double>(no_of_cycles * N_T);
    // Note: averages in do_measurement is calcuated as <g>=int d^3v fg/n, where n=int d^3v f. 

     for(size_t k{0};k<N_SPECIES;++k){
        string sp_ident = string(MATERIALS[SPECIES_KIND[k]].name);
        // calculate quantities for which ne_xt is counter_xt
        // (1) mean energy -- <epsilon>=int d^3v f*epsilon/n
        transform(meane_xt[k].begin(),meane_xt[k].end(),n_xt[k].begin(),meane_xt[k].begin(),[k](auto x, auto y){y/=FACTOR_D(k); if(y>1.0){return x/y;}else{return x;}});
        fname="meane"+sp_ident+"_xt.dat";  save_xt_1(meane_xt[k],  fname);
        // (2) mean velocity -- <u>=int d^3v f*u/n
        transform(j_xt[k].begin(),j_xt[k].end(),n_xt[k].begin(),u_xt[k].begin(),[k](auto x, auto y){y/=FACTOR_D(k); if(y>1.0){return x/y;}else{return x;}});
        fname="u"+sp_ident+"_xt.dat";  save_xt_1(u_xt[k], fname);

        // calculate P(x,t)
        transform(j_xt[k].begin(),j_xt[k].end(),efield_xt.begin(),power_xt[k].begin(),[](auto x, auto y){return x*y;});
        fname="P"+sp_ident+"_tot_xt.dat";  save_xt_1(power_xt[k], fname,FACTOR_D(k)*CHARGE(k)*SQR(INV_COUNTS_PER_BIN));

        // Do Boltzmann analysis for all species
        Boltzmann_term_analysis(k,sp_ident,INV_COUNTS_PER_BIN); 

        // calculate n(x,t) and j(x,t)
        fname="n"+sp_ident+"_xt.dat";  save_xt_1(n_xt[k], fname, INV_COUNTS_PER_BIN);
        fname="j"+sp_ident+"_xt.dat";  save_xt_1(j_xt[k], fname, FACTOR_D(k)*CHARGE(k)*INV_COUNTS_PER_BIN);
    }

    // Saving quantities in SI --> from here on, all xt-s are in SI (except for the energy)
    fname = "pot_xt.dat";     save_xt_1(pot_xt,    fname, INV_COUNTS_PER_BIN);
    fname = "efield_xt.dat";  save_xt_1(efield_xt, fname, INV_COUNTS_PER_BIN);
    // S_ion=int d^3v n_g*sigma(v)*v*f
    for (int i{0}; i<SAVE_XT_TOT; i++){
        fname = "rate_" + to_string(i) + "_xt.dat";  
        save_xt_1(process_xt[i],  fname, FACTOR_D(ELE)*INV_COUNTS_PER_BIN); 
    }
    // calculate displacement current (henceforth stored in term) as j_d=-epsilon_0*dE/dt
    differentiate(efield_xt,term,0);
    fname = "dispcurrent_xt.dat";  save_xt_1(term, fname, EPSILON0);
    Sheath_analysis();
    Current_analysis(static_cast<size_t>(N_G/2)); 
}

//---------------------------------------------------------------------//
// simulation report including stability and accuracy conditions       //
//---------------------------------------------------------------------//

void save_info(void){
    FILE     *f;
    double   plas_freq,meane,kT,debye_length,density,sim_time,e_max,v_max;
    double   coll_prob[N_SPECIES];
    double   coll_freq[N_SPECIES];
    
    density = cumul_density[ELE].at(N_G / 2) / static_cast<double>(no_of_cycles) / static_cast<double>(N_T);    // e density @ center
    plas_freq = E_CHARGE * sqrt(density / EPSILON0 / E_MASS);                                                   // e plasma frequency @ center
    meane = mean_energy_accu_center / static_cast<double>(N_center_mean_energy);                                // e mean energy @ center
    kT = 2.0 * meane * EV_TO_J / 3.0;                                                                           // k T_e @ center (approximate)
    debye_length = sqrt(EPSILON0 * kT / density) / E_CHARGE;                                                    // e Debye length @ center
    sim_time = static_cast<double>(no_of_cycles) / FREQUENCY;                                                   // simulated time
    for (int n{0}; n<N_SPECIES; n++){
        coll_freq[n] = static_cast<double>(N_coll[n]) / sim_time / static_cast<double>(x[n].size());
        coll_prob[n] = max_coll_probability[n];
    }
    
    f = fopen("info.txt","w");
    fprintf(f,"########################## PICit! simulation report ############################\n");
    fprintf(f,"Simulation parameters:\n");
    fprintf(f,"Gap distance                          = %12e [m]\n",  L);
    fprintf(f,"# of grid divisions                   = %12d\n",      N_G);
    fprintf(f,"Frequency                             = %12e [Hz]\n", FREQUENCY);
    fprintf(f,"# of time steps / period              = %12d\n",      N_T);
    fprintf(f,"# of electron / ion time steps        = %12d\n",      N_SUB);
    fprintf(f,"Voltage amplitude                     = %12e [V]\n",  VOLTAGE);
    fprintf(f,"Pressure                              = %12e [Pa]\n", PRESSURE);
    fprintf(f,"Mean temperature                      = %12e [K]\n",  mean_temp);
    fprintf(f,"perpendicular B-field                 = %12e [T]\n",  B_FIELD);
    for (int t{0}; t<N_TARGET; t++)
        fprintf(f,"initial target gas composition [%2d]   = %12.2f %% (%s)\n", t, 100.0*TARGET_RATIO[t], MATERIALS[TARGET_KIND[t]].name.data());
    for(int n{0}; n<N_SPECIES; n++)
        fprintf(f,"Superparticle weight for species %2d   = %12e (%s)\n", n, WEIGHT(n), MATERIALS[SPECIES_KIND[n]].name.data());
    fprintf(f,"# of simulation cycles in this run    = %12d\n",      no_of_cycles);
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fprintf(f,"Electron reflection coefficients      = %.2f, %.2f\n",  SURF_R_ELE[0],SURF_R_ELE[1]);
    for(int n{0}; n<N_SPECIES; n++)
        fprintf(f,"Electron emission coefficients [%2d]   = %.2f, %.2f (for %s impact) \n",  n, SURF_E_EMISSION[n][0], SURF_E_EMISSION[n][1], MATERIALS[SPECIES_KIND[n]].name.data());
    for(int n{0}; n<N_SPECIES; n++)
        fprintf(f,"Cross-section name for species %2d     = %s\n",  n, cs_names[n].data());
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fprintf(f,"Electron density @ center             = %12e [m^{-3}]\n", density);
    fprintf(f,"Plasma frequency @ center             = %12e [rad/s]\n",  plas_freq);
    fprintf(f,"Debye length @ center                 = %12e [m]\n",      debye_length);
    for (int n{0}; n<N_SPECIES; n++)
        fprintf(f,"%-12s collision frequency      = %12e [1/s]\n",   MATERIALS[SPECIES_KIND[n]].name.data(), coll_freq[n]);
    fprintf(f,"Plasma frequency @ center * DT_E      = %12.4f (OK if less than 0.20)\n", plas_freq * DT_E);
    fprintf(f,"DX / Debye length @ center            = %12.4f (OK if less than 1.00)\n", DX / debye_length);
    fprintf(f,"Electron cyclotron frequency * DT_E   = %12.4f (OK if less than 0.20)\n", B_FIELD/FACTOR_MQ(ELE) * DT_E);
    for (int n{0}; n<N_SPECIES; n++)
        fprintf(f,"Max. %-12s coll. probabability = %12.4f (OK if less than 0.05)\n", MATERIALS[SPECIES_KIND[n]].name.data(), coll_prob[n]);

    // calculate maximum electron energy for which the Courant condition holds:
    v_max = DX / DT_E;
    e_max = 0.5 * E_MASS * v_max * v_max * J_TO_EV;
    fprintf(f,"Max e- energy for CFL condition       = %12.4f [eV]\n", e_max);
    fprintf(f,"--------------------------------------------------------------------------------\n");
    for (int t{0}; t<N_TARGET; t++)
        fprintf(f,"Max. target densities: %2d             = %12e [m^{-3}]  (%s)\n", t, max_target_density[t], MATERIALS[TARGET_KIND[t]].name.data());
    fprintf(f,"--------------------------------------------------------------------------------\n");
    for(int n{0}; n<N_SPECIES; n++) {
        fprintf(f,"%-12s flux at powered electr.  = %12e [m^{-2} s^{-1}]\n", MATERIALS[SPECIES_KIND[n]].name.data(), N_arrived[n][0] * WEIGHT(n) / ELECTRODE_AREA / (no_of_cycles * PERIOD));
        fprintf(f,"%-12s flux at grounded electr. = %12e [m^{-2} s^{-1}]\n", MATERIALS[SPECIES_KIND[n]].name.data(), N_arrived[n][1] * WEIGHT(n) / ELECTRODE_AREA / (no_of_cycles * PERIOD));
        fprintf(f,"Mean %-12s energy at pow. ele. = %12e [eV]\n", MATERIALS[SPECIES_KIND[n]].name.data(), mean_energy_pow[n]);
        fprintf(f,"Mean %-12s energy at grounded  = %12e [eV]\n", MATERIALS[SPECIES_KIND[n]].name.data(), mean_energy_gnd[n]);
    }
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fprintf(f,"Flux of reflected electrons from powered electrode  = %12e [m^{-2} s^{-1}]\n", N_ele_refl[0] * WEIGHT(ELE) / ELECTRODE_AREA / (no_of_cycles * PERIOD));
    fprintf(f,"Flux of reflected electrons from grounded electrode = %12e [m^{-2} s^{-1}]\n", N_ele_refl[1] * WEIGHT(ELE) / ELECTRODE_AREA / (no_of_cycles * PERIOD));
    for(int n{0}; n<N_SPECIES; n++) {
        fprintf(f,"Flux of electrons emitted by species %2d from powered electrode  = %12e [m^{-2} s^{-1}]\n", n, N_ele_emit[n][0] * WEIGHT(ELE) / ELECTRODE_AREA / (no_of_cycles * PERIOD));
        fprintf(f,"Flux of electrons emitted by species %2d from grounded electrode = %12e [m^{-2} s^{-1}]\n", n, N_ele_emit[n][1] * WEIGHT(ELE) / ELECTRODE_AREA / (no_of_cycles * PERIOD));                
    }
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fprintf(f,"Collision rate XT files (total %d):\n", SAVE_XT_TOT);
    int XT_counter{0};
    for(int n{0}; n<N_SPECIES; n++)
        for (int i{0}; i<SAVE_XT_NUM[n]; i++){
            fprintf(f,"   rate_%d_xt.dat =  %s  (%s, %d)\n", XT_counter, sigma_name[n].at(SAVE_XT_PROCESS[n][i]).c_str(), CS_FILES[n].data(), SAVE_XT_PROCESS[n][i]);  
            XT_counter++;
        }
    fprintf(f,"--------------------------------------------------------------------------------\n");
    for(int n{0}; n<N_SPECIES; n++){
        fprintf(f,"Event numbers for %s  (%s)\n", MATERIALS[SPECIES_KIND[n]].name.data(), CS_FILES[n].data());
        for (int i{0}; i<N_CS[n]; i++){
            fprintf(f,"%6d. =  %llu  (%s)\n", i+1, event_counter.at(n).at(i), sigma_name[n].at(i).c_str());  
            XT_counter++;
        }
    }
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fclose(f);

    f = fopen("o2a_calc.txt","w");
    if (META_T_INDEX){
        fprintf(f,"## O2a density calculation ##\n");
        fprintf(f,"META_T_INDEX               = %d\n", META_T_INDEX); 
        fprintf(f,"META_SOURCE_PROCESS        = %d\n", META_SOURCE_PROCESS); 
        fprintf(f,"META_CYCLE_SKIP            = %d\n", META_CYCLE_SKIP); 
        fprintf(f,"META_DIFF_0 (1 Pa, 300K)   = %e\n", META_DIFF_0); 
        fprintf(f,"META_REFLECTION            = %e\n", META_REFLECTION); 
        double Diff_meta   = META_DIFF_0 * pow(mean_temp/300.0,1.5) / PRESSURE;
        double v_avg       = sqrt(3.0*K_BOLTZMANN*mean_temp/TARGET_MASS(META_T_INDEX));
        double lambda_meta = 2.0 * Diff_meta / v_avg;
        double beta_meta   = lambda_meta * (1.0 + META_REFLECTION) / (1.0 - META_REFLECTION) / SQRT3;
        fprintf(f,"--- derived ----\n"); 
        fprintf(f,"Diff_meta [m2/s]           = %e\n", Diff_meta); 
        fprintf(f,"v_avg [m/s]                = %e\n", v_avg); 
        fprintf(f,"lambda_meta [m]            = %e\n", lambda_meta); 
        fprintf(f,"beta_meta [m]              = %e\n", beta_meta); 
        fprintf(f,"n0_factor                  = %e\n", beta_meta / (DX + beta_meta)); 
    }
    fclose(f);
}


//------------------------------------------------------------------------------------------//
// main                                                                                     //
// command line arguments:                                                                  //
// [1]: number of cycles (0 for init)                                                       //
// [2]: "m" turns on data collection and saving                                             //
//------------------------------------------------------------------------------------------//

int main (int argc, char *argv[]){
    vector<string> args;
    int arg1;

    cout << endl << MESSAGE << endl;
    cout << ">> PICit! starting..." << endl;

    if (argc == 1) {
        cout << ">> PICit! error = need starting_cycle argument" << endl;
        return 1;
    } else {
        args.assign(argv + 1, argv + argc);
        arg1 = stoi(args.front());
        if (argc > 2) {
            if (args.at(1) == "m"){
                measurement_mode = true;                  // measurements will be done
            } else {
                measurement_mode = false;
            }
        }
    }
    if (measurement_mode) {
        cout << ">> PICit! measurement mode: on" << endl;
    } else {
        cout << ">> PICit! measurement mode: off" << endl;
    }

    read_cross_sections();
    calc_total_cross_sections();
    //for(int i{0};i<N_SPECIES;++i) { test_cross_sections(i);} return 1;

    init_vectors();

    // init move_particles function
    if (B_FIELD == 0.0) { move_particles = &move_particles_no_B; }
    else { move_particles = &move_particles_perpendicular_B; }

    

    auto start_chrono = chrono::steady_clock::now();
    
    if (arg1 == 0) {
        datafile = fopen("conv.dat","w");

 	    ofstream b("bias.dat");                           //clear bias.dat if ./PICit 0
        b<<scientific<<0.0;
        b.close();

        fill(temp.begin(), temp.end(), T_WALL);           // initialize temperature vector

        for (int i=0; i<N_TARGET; i++){                    // initialize background gas densities
            for (int k=0; k<N_G; k++){
                target_density[i].at(k) = PRESSURE / (K_BOLTZMANN * temp.at(k));
            }
            for(size_t k{1};k<N_SPECIES;++k){
                if(SPECIES_KIND[k]==TARGET_KIND[i]){
                    fill(target_density[i].begin(), target_density[i].end(),0.0);
                }
            }
        }
        calc_max_coll_probability();
        no_of_cycles = 1; 
        cycle = 1;                                        // init cycle
        init_seed(N_INIT);                                // seed initial electrons & ions
        cout << ">> PICit! running initializing cycle" << endl;
        Time = 0.0;
        old_charge     = 0.0;
        old_old_charge = 0.0;
        old_sigma      = 0.0;
        qpow           = 0.0;

        do_one_cycle();

        update_target_densities(1,0);
        cycles_done = 1;

    } else {

        datafile = fopen("conv.dat","a");

	    ifstream biasfile("bias.dat");
   	    if(!biasfile){
       	 	//if bias.dat does not exist for some reason, make one with bias=0
       		ofstream b("bias.dat");
       		b<<scientific<<0.0;
            b.close();
   	    }
    	else biasfile>>bias;
    	biasfile.close();

        no_of_cycles = arg1;                              // run number of cycles specified in command line
        load_target_densities();                          // read previous target densitites from file
        load_particle_data();                             // read previous configuration from file
        mean_temp_calc();                                 // calculate mean temperature

        cout << ">> PICit! running "<< no_of_cycles <<" cycle(s)" << endl;
        for (cycle=cycles_done+1;cycle<=cycles_done+no_of_cycles;cycle++) {
            if(SEEDING && x[ELE].size() < N_INIT){ init_seed(N_INIT,true); } 
            do_one_cycle();
            update_target_densities(cycle-cycles_done, cycle);
            double sump = 0;
            for (int k=0; k<N_G; k++){
                sump += pfield.at(k);
            }
            cout << "Average P_field: " << sump/N_G << endl;
        
        }
        cycles_done += no_of_cycles;
        solve_heat(pfield, ((double)no_of_cycles)*PERIOD, temp.front(), temp.back(), 0);
        update_target_densities(no_of_cycles, cycles_done);
        save_xvector(temp, "temp.dat"); 
        save_xvector(pfield, "pfield.dat");
        save_xvector(target_density[0], "target_density[0].dat");
    }

    fclose(datafile);
    if(alpha_0==0)
    calculate_bias(bias);
    ofstream b("bias.dat");
    b << scientific << bias; b.close();
    


    auto stop_chrono = chrono::steady_clock::now();
    
    save_particle_data();
    save_target_densities();                          
    save_density();
    if (measurement_mode) {
        cout << ">> PICit! saving diagnostics data" << endl;
        save_eepf();
        for(size_t i{0}; i<N_SPECIES; ++i){
            save_fed(i);
        }
        save_all_xt();
        save_info();
    }

    chrono::duration<double> elapsed_seconds = stop_chrono-start_chrono;
    cout << ">> PICit! simulation of " << no_of_cycles << " cycle(s) is completed in " << fixed << setprecision(2) 
         << elapsed_seconds.count() << " seconds." << endl;
}