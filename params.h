//  params_Ar.h
//  PICit! parameter file definig discharge conditions for:
//  Argon gas 

//
// Simulation Parameters
//
constexpr string_view MESSAGE = "PICit! - 1d3v electrostatic PIC/MCC plasma simulation @ Wigner RCP";

// case parameters
constexpr int         N_G             = 512;                        // number of grid points
constexpr int         N_T             = 4000;                       // time steps within an RF period
constexpr double      FREQUENCY       = 13.56e6;                     // driving frequency [Hz]
constexpr double      VOLTAGE         = 250.0;                      // voltage amplitude [V]
constexpr double      L               = 0.025;                      // electrode gap [m]
constexpr double      PRESSURE        = 10.0;                       // gas pressure [Pa]
constexpr double      T_WALL          = 300.0;                      // wall temperature
constexpr double      WEIGHT_COM      = 7.0e4;                      // common weight of superparticles (see weight factors)
constexpr double      ELECTRODE_AREA  = 1.0e-4;                     // (fictive) electrode area [m^2]
constexpr int         N_INIT          = 1000;                       // number of initial electrons and ions
constexpr int         N_SUB           = 20;                         // ions move only in these cycles (subcycling)

const int             N_EEPF  = 2000;                               // number of energy bins in Electron Energy Probability Function (EEPF)
const double          DE_EEPF = 0.05;                               // resolution of EEPF [eV]
const int             N_FED   = 200;                                // number of energy bins in Flux-Energy Distributions (FEDs)
const double          DE_FED  = 1.0;                                // resolution of FEDs [eV]

constexpr double      B_FIELD         = 0.0;                        // perpendicular magnetic field [T]

// additional (derived) constants
constexpr double      PERIOD          = 1.0 / FREQUENCY;                           // RF period length [s]
constexpr double      DT_E            = PERIOD / (double)(N_T);                    // electron time step [s]
constexpr double      DT_I            = N_SUB * DT_E;                              // ion time step [s]
constexpr double      DX              = L / (double)(N_G - 1);                     // spatial grid division [m]
constexpr double      INV_DX          = 1.0 / DX;                                  // inverse of spatial grid size [1/m]
constexpr double      GAS_DENSITY     = PRESSURE / (K_BOLTZMANN * T_WALL);         // background gas density [1/m^3]
constexpr double      OMEGA           = TWO_PI * FREQUENCY;                        // angular frequency [rad/s]
constexpr double      DV              = ELECTRODE_AREA * DX;                       // grod cell volume
constexpr double      E_WALL          = 3.0/2.0*K_BOLTZMANN*T_WALL;                // energy of wall temperatured target atom

// simulation flags
constexpr bool        SEEDING         = true;                         // additional seeding of particles to ignite plasma
constexpr bool        WARM_GAS        = false;                        // warm gas approximation for electrons --> if needed, set to true
constexpr bool        GAS_HEATING     = false;                        // gas heating calculation --> if needed, set to true
constexpr bool        COULOMB_COLL    = false;                        // Coulomb collisions --> if needed, set to true

// species identifiers
constexpr int         N_SPECIES                 = 3;                              // number of species in the simulation
constexpr int         SPECIES_KIND[N_SPECIES]   = { ELE, AR_ION, AR_FAST };
constexpr double      DT[N_SPECIES]             = { DT_E, DT_I, DT_I};          // time steps
constexpr double      WEIGHT_FACTORS[N_SPECIES] = { 1.0, 1.0, 30.0};            // superparticle weight factors
constexpr string_view CS_FILES[N_SPECIES]       = { "ar_e_cs.bin",              // cs files to use (one for each species)
                                                    "ar_arp_fa_cs.bin",
                                                    "arf_ar_cs.bin"};     

constexpr int         N_TARGET                  = 1;                             // number of collision target (background) species
constexpr int         TARGET_KIND[N_TARGET]     = { AR_GAS };                    // materials of collision target (background) species
constexpr double      TARGET_RATIO[N_TARGET]    = { 1.0 };                       // concentrations of collision target (background) species

// metastable background definition
constexpr int         META_T_INDEX        = 0;       // index of metastable target, set to 0 for skipping metastable calculation
constexpr int         META_SOURCE_PROCESS = 6;       // index of metastable source process in CS file assuming species = 0 = electrons
constexpr int         META_CYCLE_SKIP     = 20;      // No. of RF cycles between metastable density updates
constexpr double      META_DIFF_0         = 0.015;   // O2 metastable diffusion at 1 Torr, 300 K, [m2/s]
//constexpr double      META_DIFF_0         = 0.0177;  // Neon metastable diffusion at 1 Torr, 300 K, [m2/s], http://dx.doi.org/10.1088/0022-3727/49/18/185202
constexpr double      META_REFLECTION     = 0.994;

// XT measurement definition
constexpr int         SAVE_XT_NUM[N_SPECIES]         = { 0, 0, 0};               // number of processes for XT analysis (max 10 for each species)
constexpr int         SAVE_XT_PROCESS[N_SPECIES][10] = { {},{} };                // list of processes for XT analysis 
constexpr int         SAVE_XT_TOT                    = SUM_ARRAY(SAVE_XT_NUM);   // total number of process XT-s to store

// surface coefficients
constexpr int         N_SIDES              = 2;                     // No. of surfaces: 0 = powered/left, 1 = grounded/right
constexpr double      SURF_NORMAL[N_SIDES] = {1.0, -1.0};           // surface normal vector directions
constexpr double      SURF_R_ELE[N_SIDES]  = {0.7, 0.7};            // elastic electron reflection coefficients on two sides
constexpr double      SURF_E_EMISSION[N_SPECIES][N_SIDES] = { {0.0, 0.0},    // secondary electron emission yields
                                                              {0.07, 0.07},
                                                              {0.0, 0.0} };
constexpr double      ALPHA                 = 0.5;                  // Accomodation coefficient @@@ 

// Verboncoeour type solution of the Poisson equation
constexpr double      RESISTANCE   = 0.0;                           // resistor
constexpr double      INV_C        = 0.0;                           // inverse of capacitance 

// species dependent factors (for electron + heavies system)
constexpr double  CHARGE(int n)      {return MATERIALS[SPECIES_KIND[n]].charge;}     // charges
constexpr double  ABSCHARGE(int n)   {return MATERIALS[SPECIES_KIND[n]].charge < 0 ? -MATERIALS[SPECIES_KIND[n]].charge : MATERIALS[SPECIES_KIND[n]].charge;} //abscharge
constexpr double  MASS(int n)        {return MATERIALS[SPECIES_KIND[n]].mass;}       // masses
constexpr double  WEIGHT(int n)      {return WEIGHT_COM*WEIGHT_FACTORS[n];}          // superparticle weights 
constexpr double  FACTOR_D(int n)    {return WEIGHT(n) / DV;}                        // factor linking paricle number per cell to density
constexpr double  FACTOR_P(int n)    {return DT[n] / MASS(n) * CHARGE(n);}           // factor linking E-field to particle acceleration
constexpr double  FACTOR_MQ(int n)   {return MASS(n) / ABSCHARGE(n);}                // mass to charge ratios

constexpr double  TARGET_MASS(int t) {return MATERIALS[TARGET_KIND[t]].mass;}                          // target species masses
constexpr double  OPAL_W(int t)      {return MATERIALS[TARGET_KIND[t]].Opal_w;}                        // Opal w electron spectrum shape parameter
constexpr double  MU(int n, int t)   {return MASS(n)*TARGET_MASS(t) / ( MASS(n)+TARGET_MASS(t) ); }    // reduced masses

// cross sections
constexpr int         CS_RANGES            = 1000000;               // number of entries in cross section arrays
constexpr double      DE_CS                = 0.001;                 // energy division in cross section arrays [eV]
constexpr int         CS_STRING_LENGTH     = 128;                   // length of comment strings

// constants for temperature calculation
constexpr int         MAX_ITER             = 100;                        // Maximum number of iteration allowed during heat equation solving (boundary condition fitting)
const double          KAPPA                = 0.0177;                     // Thermal conductivity for argon gas [W/(m K)]
const double          cp                   = 20.7849;                    // Specific heat of Ar at const pressure [J /(K mol)]
const double          R                    = 8.31446261;                 // Universal gas constant [J /(K mol)]
const double          KHI0                 = 0.6505e-23;                 //
const double          KHI                  = (cp+1.25*R)/(cp-0.5*R);     //
const double          SIGMA_T              = 42e-20;                     // Total cross section for Ar/Ar collisions [m^2]
const double          LAMBDA               = (2-ALPHA)/ALPHA*KHI*KHI0*T_WALL/(PRESSURE*SIGMA_T);     // Constant for thermal boundary condition
const double          FA_E_THRESHOLD_FACTOR= 9.0*3.0/2.0*K_BOLTZMANN;    // Factor for fast atom thershold energy calculation