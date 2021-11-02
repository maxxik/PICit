//  consts.h
//  PICit!
//


constexpr double     PI             = 3.141592653589793;          // mathematical constant Pi
constexpr double     TWO_PI         = 2.0 * PI;                   // two times Pi
constexpr double     SQRT3          = 1.73205080756887729;        // sqrt(3)
constexpr double     SQRT2          = 1.41421356237309504;        // sqrt(2)
constexpr double     K_BOLTZMANN    = 1.38064852e-23;             // Boltzmann's constant [J/K]
constexpr double     EPSILON0       = 8.85418781e-12;             // permittivity of free space [F/m]
constexpr double     E_CHARGE       = 1.60217662e-19;             // electron charge [C]
constexpr double     E_MASS         = 9.10938356e-31;             // mass of electron [kg]
constexpr double     EV_TO_J        = E_CHARGE;                   // eV -> Joule conversion factor
constexpr double     J_TO_EV        = 1.0 / E_CHARGE;             // eV <- Joule conversion factor

// collision types
 constexpr int        COLL_ISO          = 0;                       // isotropic scattering
 constexpr int        COLL_BACK         = 1;                       // backscattering
 constexpr int        COLL_ANISO        = 2;                       // anisotropic scattering (screened-Coulomb interaction)
 constexpr int        COLL_ION          = 3;                       // ionization (e.g. Opal energy share)
 constexpr int        COLL_DISS         = 4;                       // molecular dissociation
 constexpr int        COLL_DET          = 5;                       // detachment
 constexpr int        COLL_NEUT         = 6;                       // neutralization
 constexpr int        COLL_DISS_ATTACH  = 7;                       // dissociative attachment
 constexpr int        COLL_DISS_EXC     = 8;                       // dissociative excitation
 constexpr int        COLL_DISS_REC     = 9;                       // dissociative recombination
 constexpr int        COLL_ASS_DET      = 10;                      // associative detachment
 constexpr int        COLL_ISO_FASTATOM = 11;                      // @@@
 constexpr int        COLL_BACK_FASTATOM = 12;                     // @@@
 constexpr int        COLL_CHARGE_TRANS = COLL_DISS_ATTACH;        // asymmetric charge transfer 


// particle identifiers
constexpr int        ELE             = 0;                           // electron must be always at first (0) index
constexpr int        AR_ION          = 1;
constexpr int        AR_FAST         = 2;     
constexpr int        AR_META         = 3;     
constexpr int        AR_GAS          = 4;     
constexpr int        NE_ION          = 5;     
constexpr int        NE_FAST         = 6;     
constexpr int        NE_META         = 7;     
constexpr int        NE_GAS          = 8;     
constexpr int        HE_ION          = 9;     
constexpr int        HE_FAST         = 10;     
constexpr int        HE_META         = 11;     
constexpr int        HE_GAS          = 12;     
constexpr int        O2_P_ION        = 13;
constexpr int        O2_FAST         = 14;
constexpr int        O2_META         = 15;
constexpr int        O2_GAS          = 16;     
constexpr int        O_M_ION         = 17;
constexpr int        O_GAS           = 18;
constexpr int        O3_GAS          = 19;

constexpr int        N_MATERIALS     = 20;

struct material_type{ double charge, mass, Opal_w; string_view name; };
constexpr material_type MATERIALS [N_MATERIALS] = {        // !!! order must match identifier list above !!!
    { -E_CHARGE,         E_MASS,  1.0, "electron" },          // ELE        0
    {  E_CHARGE, 6.63352090e-26, 10.0, "Argon_ion" },         // AR_ION     1
    {       0.0, 6.63352090e-26, 10.0, "Argon_fast_atom" },   // AR_FAST    2
    {       0.0, 6.63352090e-26, 10.0, "Argon_metastable" },  // AR_META    3
    {       0.0, 6.63352090e-26, 10.0, "Argon_gas" },         // AR_GAS     4
    {  E_CHARGE, 3.35091770e-26, 24.2, "Neon_ion" },          // NE_ION     5
    {       0.0, 3.35091770e-26, 24.2, "Neon_fast_atom" },    // NE_FAST    6
    {       0.0, 3.35091770e-26, 24.2, "Neon_metastable" },   // NE_META    7
    {       0.0, 3.35091770e-26, 24.2, "Neon_gas" },          // NE_GAS     8
    {  E_CHARGE, 6.64647310e-27, 15.8, "Helium_ion" },        // HE_ION     9
    {       0.0, 6.64647310e-27, 15.8, "Helium_fast_atom" },  // HE_FAST    10
    {       0.0, 6.64647310e-27, 15.8, "Helium_metastable" }, // HE_META    11
    {       0.0, 6.64647310e-27, 15.8, "Helium_gas" },        // HE_GAS     12
    {  E_CHARGE, 5.31339240e-26, 17.4, "O2p_ion" },           // O2_P_ION   13
    {       0.0, 5.31339240e-26, 17.4,  "O2_fast_molecule" },  // O2_FAST    14
    {       0.0, 5.31339240e-26, 17.4,  "O2_metastable" },     // O2_META    15
    {       0.0, 5.31339240e-26, 17.4,  "O2_gas" },            // O2_GAS     16
    { -E_CHARGE, 2.65669620e-26, 17.4,  "Om_ion" },            // O_M_ION    17
    {       0.0, 2.65669620e-26, 17.4,  "O_gas" },             // O_GAS      18
    {       0.0, 7.97008860e-26, 17.4,  "O3_gas" },            // O3_GAS     19
};


// helper constexpr definition
template <int Size> constexpr int SUM_ARRAY(const int (&arr)[Size]) {int ret = 0; for (int i=0; i<Size; ++i) ret+=arr[i]; return ret; } 
