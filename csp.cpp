// Budapest PICit! project
// cross_section_prepocessor
// 2020-2021

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <array>
#include <limits>

using namespace std;

// constants (shared with the simulation)
#include "consts.h"

// Cross-section preprocesor parameters (shared with the simulation)
#include "params.h"

const string CS_FOLDER = "cross_sections";

struct cross_section_data_type {
  int    type;                    // collision type
  int    target;                  // target particle identifier
  int    product[4];              // list of product particles
  double E_threshold;             // threshold energy
  std::string process;            // process name string
  std::vector<float> values;      // actual cross-section values in SI units
};


//
//  Global variables
//
int      cs_sets;      // # of CS sets (processes) included

string   interpolation_type;
string   input_name;
string   output_name;
string   comments;

vector<cross_section_data_type> sigma;

typedef double  cross_section [CS_RANGES+1];
cross_section   detach1, detach2;

int target_particle;


//
// reads LxCat format cross section data file
//
void read_cs_LxCat(const string filename, int& sets, bool notonlyread=true, int cs_nr=0, double reduced_mass_factor=1.0){
  ifstream       file;
  string         line;
  vector<string> lines;
  vector<int>    dash_lines;
  vector<double> en_v, cs_v;
  int            i, j, k, start_line, stop_line, number_of_lines;
  double         en, cs_value;
  
  // read in all lines from file
  
  cout << "Reading cross section --> " << filename << endl;
  file.open(filename);
  if (!file.is_open()) {                                              // stop on error
    cerr << "ERROR_1: file not found"  << endl; throw exception(); }
  
  while (getline(file,line)){                                         // read file line-by-line
    lines.push_back(line);                                            // add each line to "lines"
  }
  file.close();
  number_of_lines = lines.size();
  
  // find lines with dashes at beginning, list their indexes in "dash_lines"
  
  if(notonlyread){
  for (j=0; j<int(lines.size()); j++) {
    if (!lines.at(j).empty() && lines.at(j)[lines.at(j).size() - 1] == '\r')  // remove last \r character
      lines.at(j).erase(lines.at(j).size() - 1);
    if (lines.at(j).rfind("-----", 0) == 0)               // searching only at the beginning of the line
      dash_lines.push_back(j);
  }
  
  // this many c.s. sets have been found, allocate memory:
  
  sets = dash_lines.size() / 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
  }
  
  // find types of processes, process names and values
    
  // ZD: for ELASTIC collisions LxCat gives the mass ratio in the same line as
  // the threshold energy for inelastic processes. We do not use this number
  // and put zero to the E_threshold field, which is true for an elastic process.
  // The mass ratio can be taken from the declaration of the parameters of species.
    
  k=0;
  for (j=0; j<number_of_lines; j++) {
    if (lines.at(j).rfind("ELASTIC", 0) == 0) {
      sigma.at(k).type = COLL_ISO;
      sigma.at(k).process = lines.at(j+1);
      sigma.at(k).E_threshold = 0;                 // setting of zero energy here
      k++;
    }
    if (lines.at(j).rfind("BACKSCATTERING", 0) == 0) {
      sigma.at(k).type = COLL_BACK;
      sigma.at(k).process = lines.at(j+1);
      sigma.at(k).E_threshold = 0;                 // setting of zero energy here
      k++;
    }
    if (lines.at(j).rfind("EXCITATION", 0) == 0) {
      sigma.at(k).type = COLL_ISO;
      sigma.at(k).process = lines.at(j+1);
      sigma.at(k).E_threshold = stof(lines.at(j+2));
      k++;
    }
    if (lines.at(j).rfind("IONIZATION", 0) == 0) {
      sigma.at(k).type = COLL_ION;
      sigma.at(k).process = lines.at(j+1);
      sigma.at(k).E_threshold = stof(lines.at(j+2));
      k++;
    }
  }
  
  if (k != sets){ cout << "Warning: inconsistent No. of data-sets found!" << endl; throw exception(); }  // stop
  }


  double a, b;
  int st;
  if(!notonlyread) st=1; else{st=sets;}
  for (k=0; k<st; k++){
    
    // read in values for given set and store in en_v and cs_v
    
    en_v.clear();
    cs_v.clear();
    if(notonlyread){
      start_line = dash_lines.at(2 * k);
      stop_line  = dash_lines.at(2 * k + 1);
      for (j=start_line+1; j<=stop_line-1; j++){
        std::istringstream iss(lines.at(j));
        if (!(iss >> a >> b)) { cout << "Conversion error" << endl; throw exception(); }   // stop on error
        en_v.push_back(a);
        cs_v.push_back(fmax(b,std::numeric_limits<float>::min()));     /////////////
      }
    }
    else{
      start_line = 0;
      stop_line  = number_of_lines;
      for (j=start_line; j<stop_line; j++){
        std::istringstream iss(lines.at(j));
        if (!(iss >> a >> b)) { cout << "Conversion error" << endl; throw exception(); }   // stop on error
        en_v.push_back(a);
        cs_v.push_back(fmax(b,std::numeric_limits<float>::min()));     /////////////
      }
    }

    
    
    // LOG or LIN interpolation:
    
    if (interpolation_type.compare("LOG") == 0) {          // LOG
      for(i=0; i<CS_RANGES; i++){
        en = DE_CS * max(i, 1) * reduced_mass_factor;
        if (en <= en_v.front()){
           if(notonlyread){cs_value = cs_v.front();} else{cs_value=std::numeric_limits<float>::min();}
           }
        else if (en >= en_v.back()) cs_value  = cs_v.back();
        else {
             j=1;
            while (en >= en_v.at(j)) { j++; }
            j--;
            if(en_v.at(j)==0)en_v.at(j)=1e-10;
	          if(cs_v.at(j)==0)cs_v.at(j)=std::numeric_limits<float>::min();;

            cs_value = cs_v.at(j)*exp( log(cs_v.at(j+1)/cs_v.at(j)) * log(en/en_v.at(j)) / log(en_v.at(j+1)/en_v.at(j)) );
          }
        
        if(notonlyread){ sigma.at(k).values.at(i) = (float) cs_value;}
        else{//overwrite CS
          sigma.at(cs_nr).values.at(i) = (float) cs_value;}
      }
    }
    else {                                                 // LIN
      for(i=0; i<CS_RANGES; i++){
        en = DE_CS * max(i, 1) * reduced_mass_factor;
        if (en <= en_v.front()){
           if(notonlyread){cs_value = cs_v.front();} else{cs_value=std::numeric_limits<float>::min();}
           }
	      else if(en >= en_v.back()){ cs_value  = cs_v.back();} 
	      else {
            j=1;
            while (en >= en_v.at(j)) { j++; }
            j--;
            cs_value = cs_v.at(j) + (cs_v.at(j+1) - cs_v.at(j)) * (en - en_v.at(j)) / (en_v.at(j+1) - en_v.at(j));
        }
        if(notonlyread){ sigma.at(k).values.at(i) = (float) cs_value;}
        else{//overwrite cs
          sigma.at(cs_nr).values.at(i) = (float) cs_value;}
       
      }
    }
  }

  if(notonlyread){
  (void)cs_nr;
    for (k=0; k<sets; k++){
      cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
    }
  }
}


//
// Merge two binary CS files into one
//
void merge_CS_files(string dest, string src1, string src2){
  vector<char>  buffer1;
  vector<char>  buffer2;
  char          st [CS_STRING_LENGTH];
  ifstream      file;
  ofstream      ofile;
  double        cs_de, ds_de;
  int           cs_n, cs_sets, ds_n, ds_sets, b1size, b2size;
  string        csname, interpolation_type, comment1, comment2;

  csname.reserve(CS_STRING_LENGTH);
  interpolation_type.reserve(CS_STRING_LENGTH);
  comment1.reserve(CS_STRING_LENGTH);
  comment2.reserve(CS_STRING_LENGTH);
  cout << ">> CS preprocessor: Merging cross-sections files " << src1 << " and " << src2 << " into " << dest << endl;

  // loading and printing src1 information
  cout << "   Loading cross-sections from file " << src1 << endl;
  file.exceptions(ifstream::failbit | ifstream::badbit);
  try { file.open(src1, ios::binary); }
  catch (system_error& e) { cerr << e.code().message() << endl; }
  file.read(st, CS_STRING_LENGTH); csname = st;      
  file.read(st, CS_STRING_LENGTH); interpolation_type = st;
  file.read(st, CS_STRING_LENGTH); comment1 = st;
  file.read((char*) &cs_de,   sizeof(double));
  file.read((char*) &cs_n,    sizeof(int));
  file.read((char*) &cs_sets, sizeof(int));
  if(cs_n  != CS_RANGES) cout << "Warning: Cross-section data lentgh inconsistent !!!" << endl;
  if(cs_de != DE_CS)     cout << "Warning: Cross-section data resolution inconsistent !!!" << endl;
  cout << "      Cross-section name:  " << csname << endl;
  cout << "      Interpolation type:  " << interpolation_type << endl;
  cout << "      Additional comments: " << comment1 << endl;
  cout << "      Energy resolution:   " << cs_de << " eV" << endl;
  cout << "      Total data sets:     " << cs_sets << endl;
  cout << "      Data points per set: " << cs_n << endl;
  b1size = cs_sets * (cs_n*sizeof(float) + sizeof(double) + 6*sizeof(int) + CS_STRING_LENGTH);
  buffer1.resize(b1size);
  file.read((char*) buffer1.data(), b1size);
  file.close();

  // loading and printing src2 information
  cout << "   Loading cross-sections from file " << src2 << endl;
  file.exceptions(ifstream::failbit | ifstream::badbit);
  try { file.open(src2, ios::binary); }
  catch (system_error& e) { cerr << e.code().message() << endl; }
  file.read(st, CS_STRING_LENGTH);       
  cout << "      Cross-section name:  " << st << endl;
  file.read(st, CS_STRING_LENGTH); 
  cout << "      Interpolation type:  " << st << endl;
  file.read(st, CS_STRING_LENGTH); comment2 = st;
  cout << "      Additional comments: " << st << endl;
  file.read((char*) &ds_de,   sizeof(double));
  file.read((char*) &ds_n,    sizeof(int));
  file.read((char*) &ds_sets, sizeof(int));
  if(ds_n  != CS_RANGES) cout << "Warning: Cross-section data lentgh inconsistent !!!" << endl;
  if(ds_de != DE_CS)     cout << "Warning: Cross-section data resolution inconsistent !!!" << endl;
  cout << "      Energy resolution:   " << ds_de << " eV" << endl;
  cout << "      Total data sets:     " << ds_sets << endl;
  cout << "      Data points per set: " << ds_n << endl;
  b2size = ds_sets * (ds_n*sizeof(float) + sizeof(double) + 6*sizeof(int) + CS_STRING_LENGTH);
  buffer2.resize(b2size);
  file.read((char*) buffer2.data(), b2size);
  file.close();

  // saving combined CS file
  cout << "   Saving cross-sections to file " << dest << endl;
  cs_sets += ds_sets;
  comment1.append(" AND ");
  comment1.append(comment2);
  comment1.resize(CS_STRING_LENGTH);
  ofile.open (dest, ios::binary);
  csname = "combined pre-processed cross section data";
  ofile.write((char*) csname.data(), CS_STRING_LENGTH);
  ofile.write((char*) interpolation_type.data(), CS_STRING_LENGTH);
  ofile.write((char*) comment1.data(), CS_STRING_LENGTH);
  ofile.write((char*) &cs_de,          sizeof(double) );
  ofile.write((char*) &cs_n,           sizeof(int) );
  ofile.write((char*) &cs_sets,        sizeof(int) );
  ofile.write((char*) buffer1.data(), b1size );
  ofile.write((char*) buffer2.data(), b2size );
  ofile.close();
}


//
// fill data with Phelps formulas for ele + Ar
//
void fill_cs_Phelps_ele(int& sets){
  
  int    k;
  double en, cs_value;
  sets = 3;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
  }

  // elastic (isotropic)
  k=0;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).process = "e + Ar -> e + Ar  isotropic elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    cs_value = (fabs(6.0/exp(3.3*log(1.0+(en/0.1)+(en/0.6)*(en/0.6)))-
               1.1*exp(1.4*log(en))/(1.0 + exp(1.2*log(en/15.0)))/
               sqrt(1+exp(2.5*log(en/5.5)) + exp(4.1*log(en/60.0))))+0.05/(1.0+en/10.0)/(1.0+en/10.0)
               + 0.01*en*en*en/(1.0+exp(6.0*log(en/12.0)))) * 1e-20;
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  
  // excitation (isotropic)
  k=1;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).process = "e + Ar -> e + Ar*  isotropic excitation";
  sigma.at(k).E_threshold = 11.5;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    if (en > 11.5){
      cs_value = 0.85*(4e-22*exp(1.1*log(en-11.5))*(1+exp(2.8*log(en/15.0)))/
              (1.0+exp(5.5*log(en/23.0)))
              +2.7e-22*(en-11.5)/exp(1.9*log(1+(en/80.0))));
    } else {
      cs_value = 0;
    }
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  
  // ionization (isotropic)
  k=2;
  sigma.at(k).type = COLL_ION;
  sigma.at(k).process = "e + Ar -> e + e + Ar+  ionization";
  sigma.at(k).E_threshold = 15.8;
  sigma.at(k).product[0] = ELE;
  sigma.at(k).product[1] = AR_ION;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    if (en>15.8){
      cs_value = 9.7e-18*(en-15.8)/(70+en)/(70+en)+6e-22*(en-15.8)*(en-15.8)*exp(-en/9);
    } else {
      cs_value = 0;
    }
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  
  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }
}

// ZD: for ions we always store the cross sections as a function of the com energy
// whenever these data are originally given as a function of the lab energy, the
// conversion from lab to com energy is to be carried out in this preprocessor!


// ZD: fill data with Phelps formulas for He+ + He from jila.colorado.edu
// Unfortunately JILA does not provide anymore these data.
// These formulas are given in the com frame, straightforward to convert

void fill_cs_Phelps_website_Hep(int& sets){
  
  const double REDUCED_MASS_RATIO = 2.0;
  int    k;
  double en, cs_value;
  sets = 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = HE_GAS;
}

  // elastic (isotropic)
  k=0;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).process = "He+ + He -> He+ + He  isotropic elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab
    cs_value = 7.63 / sqrt(en) * 1e-20;
    sigma.at(k).values.at(i) = (float) cs_value;
  }

  // elastic (backward) = symmetric charge transfer
  k=1;
  sigma.at(k).type = COLL_BACK;
  sigma.at(k).process = "He+ + He -> He+ + He  backward elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab
    cs_value = 1.0e-19 * pow(en*0.001,-0.15) * pow(1.0+0.001*en,-0.25) * pow(1.0+5.0/en,-0.15);
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  
  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }
}



// ZD: fill data with Phelps formulas for Ar+ + Ar from jila.colorado.edu
// Unfortunately JILA does not provide anymore these data.
// These formulas are given in the com frame, straightforward to convert

void fill_cs_Phelps_website_Arp(int& sets){
  
  const double REDUCED_MASS_RATIO = 2.0;
  int    k;
  double en, cs_value;
  sets = 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
  }

  // elastic (isotropic)
  k=0;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).process = "Ar+ + Ar -> Ar+ + Ar  isotropic elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab
    cs_value = 21.6 / sqrt(en) * 1e-20;
    sigma.at(k).values.at(i) = (float) cs_value;
  }

  // elastic (backward) = symmetric charge transfer
  k=1;
  sigma.at(k).type = COLL_BACK;
  sigma.at(k).process = "Ar+ + Ar -> Ar+ + Ar  backward elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab
    cs_value = 52.0 / pow(en,0.08) / (1+0.08/en) / pow(1.0+en/1000,0.3) * 1e-20;
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  
  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }
}


// ZD: fill data with Phelps formulas for Ar+ + Ar from J. Appl. Phys. 76, 747 (1994)
// These formulas are given in the lab frame, energy values need to be converted

void fill_cs_Phelps_1994paper_Arp(int& sets){
  
  const double REDUCED_MASS_RATIO = 2.0;
  int    k;
  double en, cs_value;
  sets = 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
  }

  // elastic (isotropic)
  k=0;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).process = "Ar+ + Ar -> Ar+ + Ar  isotropic elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab
    cs_value = 2e-19 * pow(en,-0.5) / (1.0+en) + 3e-19 * en / pow(1.0 + en/3.0, 2);
    sigma.at(k).values.at(i) = (float) cs_value;
  }

  // elastic (backward) = symmetric charge transfer
  k=1;
  sigma.at(k).type = COLL_BACK;
  sigma.at(k).process = "Ar+ + Ar -> Ar+ + Ar  backward elastic";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab
    cs_value = (1.15e-18 * pow(en,-0.1) * pow(1.0+0.015/en, 0.6) -
               (2e-19 * pow(en,-0.5) / (1.0+en) + 3e-19 * en / pow(1.0 + en/3.0, 2))) / 2.0;      // (Qm - Qi)/2
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  
  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }
}


//
// Writing log file for documentation
//

void print_log(void){
  FILE *f;
  int  k;
  static int number = 0;
  
  std::string name = "PICit_CSP_"+ to_string(number)+".log";
  f = fopen(name.c_str(), "w");
  fprintf(f,"--------------------------------------------------\n");
  fprintf(f,"--- PICit! Cross Section Preprocessor log file ---\n");
  fprintf(f,"--------------------------------------------------\n");
  fprintf(f,"Input file         : %s\n",   input_name.data() );
  fprintf(f,"Output file        : %s\n",   output_name.data() );
  fprintf(f,"Comments           : %s\n",   comments.data() );
  fprintf(f,"Energy division    : %6.3f\n",DE_CS);
  fprintf(f,"Energy entries     : %6d\n",  CS_RANGES);
  fprintf(f,"Interpolation      : %s\n",   interpolation_type.data() );
  fprintf(f,"Reactions          : \n");
  for (k=0; k<cs_sets; k++){
    fprintf(f,"%d: >  %s / type : %d / target : %d / value : %12e\n",k,sigma.at(k).process.data(), sigma.at(k).type, sigma.at(k).target, sigma.at(k).E_threshold);
  }
  fprintf(f,"--------------------------------------------------\n");
  fclose(f);
  ++number;
}


//
// writing to file
//
void write_cs(const string& filename, const string& comments){
  std::ofstream  file;
  std::string    csname;
  double         d = DE_CS;
  int            n = CS_RANGES;

  csname.reserve(CS_STRING_LENGTH);
  
  file.open (filename, ios::binary);
  csname = "pre-processed cross section data";
  file.write((char*) csname.data(), CS_STRING_LENGTH);
  file.write((char*) interpolation_type.data(), CS_STRING_LENGTH);
  file.write((char*) comments.data(), CS_STRING_LENGTH);
  file.write((char*) &d,             sizeof(double) );
  file.write((char*) &n,             sizeof(int) );
  file.write((char*) &cs_sets,       sizeof(int) );
  for (int k = 0; k < cs_sets; k++){
    file.write((char*) &sigma.at(k).type,            sizeof(sigma.at(k).type) );
    file.write((char*) &sigma.at(k).target,          sizeof(sigma.at(k).target) );
    file.write((char*) &sigma.at(k).E_threshold,     sizeof(sigma.at(k).E_threshold) );
    file.write((char*) sigma.at(k).product,          4*sizeof(int) );
    file.write((char*) sigma.at(k).process.data(),   CS_STRING_LENGTH );
    file.write((char*) sigma.at(k).values.data(),    CS_RANGES*sizeof(float) );
  }
  file.close();

  print_log();
}


//
// Entry point for ELE + He from Biagi v7.1 @ LxCat
//
void process_ele_He_Biagi(void){
  interpolation_type = "LIN";
  comments           = "ELE + He complete set from Biagi v7.1 @ LxCat";
  input_name         = CS_FOLDER + "/Biagi71-helium.txt";
  output_name        = "he_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = HE_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = HE_ION; };
  }

}

//
// Entry point for ELE + Ar from Biagi v7.1 @ LxCat
//
void process_ele_Ar_Biagi(void){
  interpolation_type = "LOG";
  comments           = "ELE + Ar complete set from Biagi v7.1 @ LxCat";
  input_name         = CS_FOLDER + "/Biagi71-argon.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for ELE + Ar from BSR-500 @ LxCat
//
void process_ele_Ar_BSR(void){
  interpolation_type = "LOG";
  comments           = "ELE + Ar complete set from BSR-500 @ LxCat";
  input_name         = CS_FOLDER + "/BSR-500_2014_argon.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for ELE + Ar from COP @ LxCat
//
void process_ele_Ar_COP(void){
  interpolation_type = "LOG";
  comments           = "ELE + Ar complete set from COP (Complex Optical Potential) @ LxCat";
  input_name         = CS_FOLDER + "/COP-argon.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for ELE + Ar from Hayashi @ LxCat
//
void process_ele_Ar_Hayashi(void){
  interpolation_type = "LOG";
  comments           = "ELE + Ar complete set from Hayashi @ LxCat";
  input_name         = CS_FOLDER + "/Hayashi-argon.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for ELE + Ar from Morgan (Kinema Research & Software) @ LxCat
//
void process_ele_Ar_Morgan(void){
  interpolation_type = "LOG";
  comments           = "ELE + Ar complete set from Morgan @ LxCat";
  input_name         = CS_FOLDER + "/Morgan-argon.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for ELE + Ar from Puech @ LxCat
//
void process_ele_Ar_Puech(void){
  interpolation_type = "LOG";
  comments           = "ELE + Ar complete set from Puech @ LxCat";
  input_name         = CS_FOLDER + "/Puech-argon.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for He+ + He from Phelp database (hosted by LxCat)
//
void process_He_He_Phelps(void){
  interpolation_type = "LIN";
  comments           = "He+ + He elastic set from Phelp formulas (hosted by LxCat)";
  input_name         = "";
  output_name        = "he_hep_cs.bin";
  fill_cs_Phelps_website_Hep(cs_sets);
}

//
// Entry point for ELE + AR from LISBON @ LxCat
//
void process_ele_Ar_LISBON(void){
  interpolation_type = "LOG";
  comments           = "ELE + AR complete set from LISBON @ LxCat";
  input_name         = CS_FOLDER + "/LISBON_AR_CS.txt";
  output_name        = "ar_e_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = AR_ION; };
  }
}

//
// Entry point for ELE + AR from Phelps formulas
//
void process_ele_Ar_Phelps(void){
  interpolation_type = "FUNCTION";
  comments           = "ELE + AR complete set based on Phelps formulas";
  input_name         = "";
  output_name        = "ar_e_cs.bin";
  fill_cs_Phelps_ele(cs_sets);
}

//
// Entry point for AR+ + AR from Phelps formulas from JILA website
//
void process_Arp_Ar_Phelps_website(void){
  interpolation_type = "FUNCTION";
  comments           = "AR+ + AR complete set based on Phelps formulas (JILA website)";
  input_name         = "";
  output_name        = "ar_arp_cs.bin";
  fill_cs_Phelps_website_Arp(cs_sets);
}

//
// Entry point for AR+ + AR from Phelps formulas from 1994 J. Appl. Phys. paper
//
void process_Arp_Ar_Phelps_1994paper(void){
  interpolation_type = "FUNCTION";
  comments           = "AR+ + AR complete set based on Phelps formulas (1994 JAP paper)";
  input_name         = "";
  output_name        = "ar_arp_cs.bin";
  fill_cs_Phelps_1994paper_Arp(cs_sets);
}

//
// Entry point for AR fast + AR fast from @@@ TODO
//
void process_Arf_Ar(void){
  interpolation_type = "LOG";
  comments           = "Ar fast + AR complete set from TODO @ LxCat";
  input_name         = CS_FOLDER + "/fastAr_Ar_CS.txt";
  output_name        = "arf_ar_cs.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = AR_GAS;
  }
}


//
// Entry point for ELE excit Ne2p1 state from Biagi v7.1 @ LxCat
//
void process_ele_Ne2p1_excit_Biagi(void){
  interpolation_type = "LIN";
  comments           = "ELE + Ne excit. to Ne2p1 state from Biagi v7.1 @ LxCat";
  input_name         = CS_FOLDER + "/Ne2p1excit_LXCAT_Biagi-v7.1.txt";
  output_name        = "cs_e_excit_Ne2p1_Biagi-v7.1.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = NE_META; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = NE_GAS;
  }
}

//
// Entry point for ELE excit Ne2p1 state from Puech @ LxCat
//
void process_ele_Ne2p1_excit_Puech(void){
  interpolation_type = "LIN";
  comments           = "ELE + Ne excit. to Ne2p1 state from Puech @ LxCat";
  input_name         = CS_FOLDER + "/Ne2p1excit_LXCAT_Puech.txt";
  output_name        = "cs_e_excit_Ne2p1_Puech.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = NE_META; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = NE_GAS;
  }
}

//
// Entry point for ELE excit Ne2p1 state from BSR @ LxCat
//
void process_ele_Ne2p1_excit_BSR(void){
  interpolation_type = "LIN";
  comments           = "ELE + Ne excit. to Ne2p1 state from BSR @ LxCat";
  input_name         = CS_FOLDER + "/Ne2p1excit_LXCAT_BSR.txt";
  output_name        = "cs_e_excit_Ne2p1_BSR.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = NE_META; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = NE_GAS;
  }
}

//
// Entry point for ELE + Ne (Biagi v7.1) @ LxCat
//
void process_ele_Ne_Biagi71(void){
  interpolation_type = "LIN";
  comments           = "ELE + Ne (Biagi v7.1) @ LxCat";
  input_name         = CS_FOLDER + "/Biagi71-neon.txt";
  output_name        = "cs_Ne_e.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = NE_GAS;
    if (sigma.at(k).type == COLL_ION){ sigma.at(k).product[0] = ELE; sigma.at(k).product[1] = NE_ION; };
  }
}

//
// Entry point for Ne+ + Ne (ISO+CHT) Ref: Jovanovic, Vrhovac, Petrovic, Eur. Phys. J. D 21, 335-342 (2002)
//
void process_Nep_Ne_Jovanovic(void){
  interpolation_type = "LIN";
  comments           = "Ne+ + Ne [Jovanovic, et.al. EPJ-D 21, 335-342 (2002)]";
  input_name         = CS_FOLDER + "/Ne_Ne_Jovanovic.txt";
  output_name        = "cs_Ne_Nep.bin";
  read_cs_LxCat(input_name, cs_sets);
  for (int k=0; k<cs_sets; k++){
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
    sigma.at(k).target = NE_GAS;
  }
}



//
// fill data with xpdp1 cross sections for ele + O2
//  
void fill_cs_O2_ele(int& sets){
  
  int    k;
  double en, cs_value, sigma1;
  double temp1, temp2, temp3, temp4, temp5, temp6;

  comments           = "CS set for ELE + O2 / O- / O2+ ";
  output_name        = "cs_O2_e.bin";

  sets = 17;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
  }
  

  // elastic 
  k=0;
  sigma.at(k).type =  COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2 -> e + O2             elastic scattering";
  sigma.at(k).E_threshold = 0.0;
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    if(0.0 <= en && en <= 1.2){
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      
      temp5 = -12.6132*temp3 +39.2258*temp2 -43.3875*temp1 +23.5722*en +.4464;
      temp5 *= 1e-20;
    } else if (1.2 < en && en <= 20.0){
        temp1 = en*en;
        temp2 = en*temp1;
        temp3 = en*temp2;
        temp4 = en*temp3;
      
        temp5 = -4.0554e-5*temp4 +2.7604e-3*temp3 -.07107*temp2 +.82961*temp1 -3.9163*en +11.735;
        temp5 *= 1e-20;
    } else if(20.0 < en && en <= 100.0){
      temp1 = en*en;
      temp5 = 1.3874e-4*temp1 -.0417*en +9.254364;
      temp5 *= 1e-20;
    } else temp5 = 6.5e-20;

    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;
  }

  // overwrite data
  // use elastic momentum transfer cross section and isotropic scattering
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/e_elastic_scattering_O2-Biagi.dat",cs_sets,false,0); 


  // rotational excitation (E loss = 0.02 eV)          
  k=1;
  sigma.at(k).type =  COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2(r=0) -> e + O2(r>0)   rotational excitation";
  sigma.at(k).E_threshold = 0.02;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp6 = 0.0;
    if(0.07 <= en && en <= 1.67){
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
      temp5 = en*temp4;
      temp6 = -.0859*temp5 +.4233*temp4 -.7366*temp3 +.5205*temp2 -.1537*temp1 +.0604*en -.0022;
      temp6 *= 1e-20;
    }
    cs_value = temp6;  
    sigma.at(k).values.at(i) = (float) cs_value;
  }

  // vibrational excitation v=1 (E loss = 0.19 eV)         
  k=2;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2(v=0) -> e + O2(v=1)   vibrational excitation";
  sigma.at(k).E_threshold = 0.19; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp6= 0.0;
    if(0.19 <= en && en <= 1.0) {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
    
      temp6 = -1.508395*temp4 +6.521786*temp3 -9.574636*temp2 +5.092031*temp1 -0.41602*en -0.066398;
      temp6 *= 1e-20;
      if(temp6 <= 0.0) temp6=0.0;
    } 
    else if (1.0 < en && en <= 1.67) temp6 = -7.2193e-22*(en -1.67);
    else if(4.0 <= en && en <= 5.5) temp6 = 4.572852e-22*(en -4.0);
    else if(5.5 < en && en <= 16.0){
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
      temp5 = en*temp4;
    
      temp6 = 1.0835e-6*temp5 -9.229e-5*temp4 +0.0030853*temp3 -0.050981*temp2 +0.427934*temp1 -1.6682*en +2.3919;
      temp6 *= 1e-20;
      if(temp6 <= 0.0) temp6=0.0;
    }
    else if(16.0 < en && en <= 25.0) temp6 = -4.098144e-23*(en -25.0);
  
    cs_value = temp6;  
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // vibrational excitation v=2 (E loss = 0.38 eV)         
  k=3;
  sigma.at(k).type =  COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2(v=0) -> e + O2(v=2)   vibrational excitation";
  sigma.at(k).E_threshold = 0.38; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5 = 0.0;
    if(.38 <= en && en <= 1.67)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      
      temp5 = -0.606022*temp3 +3.157773*temp2 -5.933895*temp1 +4.664064*en -1.233443;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5 = 0.0;
    }
    else if(4.0 <= en && en <= 14.5)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
      
      temp5 = -3.1339358e-6*temp4 +2.0994236e-4*temp3 -.00503577*temp2 +0.0515*temp1 -0.2074798*en +0.279;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(14.5 < en && en <= 25.0) temp5 = -1.71326e-23*(en -25.0);
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // vibrational excitation v=3 (E loss = 0.57 eV)      
  k=4;
  sigma.at(k).type =  COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2(v=0) -> e + O2(v=3)   vibrational excitation";
  sigma.at(k).E_threshold = 0.57; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5 = 0.0;
    if(.57 <= en && en <= 1.67)
    {
      temp1 = en*en;
      
      temp5 = -.055083*temp1 +.12457*en -.057531;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(4.0 <= en && en <= 15.0)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
      
      temp5 = -7.969385e-6*temp4 +4.78119632e-4*temp3 -.0107124*temp2 +0.1095564*temp1 -0.4962553*en +0.80444;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(15.0 < en && en <= 20.0) temp5 = -1.76e-23*(en -20.0);

    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  
  // vibrational excitation v=4  (E loss = 0.75 eV)      
  k=5;
  sigma.at(k).type =  COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2(v=0) -> e + O2(v=4)   vibrational excitation";
  sigma.at(k).E_threshold = 0.75;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp6= 0.0;
    if(.75 <= en && en <= 0.85) temp6 = 2.795e-25*(en - 0.75);
    else if(.85 <= en && en <= 1.67)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      
      temp6 = -0.049346*temp2 +0.16616*temp1 -0.174061*en +0.058213;
      temp6 *= 1e-20;
      if(temp6 <= 0.0) temp6=0.0;
    }
    else if(6.0 <= en && en <= 15.0)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      
      temp6 = -1.3846154e-5*temp3 +8.8449e-4*temp2 -0.020271*temp1 +0.19111*en -0.589505;
      temp6 *= 1e-20;
      if(temp6 <= 0.0) temp6=0.0;
    }

    cs_value = temp6;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // O2 SINGLET DELTA metastable excitation  (E LOSS = 0.977 eV)     
  k=6;
  sigma.at(k).type =  COLL_ISO;   
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = O2_META;
  sigma.at(k).process = "e + O2 -> e + O2(a1Delta)    metastable excitation (0.977 eV)";
  sigma.at(k).E_threshold = 0.977;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5 = 0.0;
    if(0.977 <= en && en <= 10.5)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
      
      temp5 = -3.0913e-6*temp4 +1.436e-4*temp3 -.0022876*temp2 +.0133286*temp1 -.0100266*en -.0015636;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(10.5 < en && en <= 45.0)
    {
      temp1 = en*en;
      temp2 = en*temp1;

      temp5 = -1.0959e-6*temp2 +1.349e-4*temp1 -.005984*en +.1079;
      temp5 *= 1e-20;
    }
    else if(45 < en && en <= 100.0)
    {
      temp5 = -2.18e-4*en +2.18e-2;
      temp5 *= 1e-20;
    }

    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  
  //  O2 B SINGLET SIGMA metastable excitation (isotropic) (E LOSS = 1.627 eV)    
  k=7;
  sigma.at(k).type =  COLL_ISO;  
  sigma.at(k).target = O2_GAS; 
  sigma.at(k).product[0] = O2_META;
  sigma.at(k).process = "e + O2 -> e + O2(b1Sigma)    metastable excitation (1.63 eV)";
  sigma.at(k).E_threshold = 1.627;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5 = 0.0;
    if(1.627 <= en && en <= 25.0)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;
      temp4 = en*temp3;
      
      temp5 = 5.34245e-8*temp4 -4.7117e-6*temp3 +1.581e-4*temp2 -2.4783e-3*temp1 +1.70373e-2*en -2.2343e-2;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(25.0 < en && en <= 45.0)
    {
      temp1 = en*en;
      temp2 = en*temp1;

      temp5 = 5.445e-7*temp2 -5.674e-5*temp1 +1.7502e-3*en -1.0375e-2;
      temp5 *= 1e-20;
    }
    else if(45 < en && en <= 100.0)
    {
      temp5 = -5.64e-5*en +5.64e-3;
      temp5 *= 1e-20;
    }
  
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  //  dissociative attachment  (E_t=4.2 eV)    
  k=8;
  sigma.at(k).type = COLL_DISS_ATTACH;  
  sigma.at(k).target = O2_GAS; 
  sigma.at(k).product[0] = O_M_ION;
  sigma.at(k).product[1] = O_GAS;
  sigma.at(k).process = "e + O2 -> O + O-             dissociative attachment (E_t=4.2 eV)";
  sigma.at(k).E_threshold = 4.2;   
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    if(en <= 7.0)
     detach1[i] = 1.4064e-22*exp(-(en -6.5)*(en -6.5)/1.1766);
    else if(7.0 <en && en <=15.0)
      detach1[i] = 9.0e-24 -6.0e-25*en +1.4064e-22*exp(-(en -6.5)*(en -6.5)/1.1766);
    if(en <15.0)
      detach2[i] = 0.0;
    else if(15.0<=en && en <31.0) 
      detach2[i] = 4.66e-23*exp(-(en -30.0)*(en -30.0)/50.0);
    else
      detach2[i] = 5.1745e-23 -1.96e-25*en; 
    if(en < 15.0) {
      sigma1 = detach1[i];
    }
    else {
      sigma1 = detach2[i];
    }

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // excitation  (E LOSS = 4.5 eV)    
  k=9;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "e + O2 -> e + O2             excitation (4.5 eV)";
  sigma.at(k).E_threshold = 4.5;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5= 0.0;
    if(4.5 <= en && en <= 4.8)
    {
      temp5 = 0.01*en -.045;
      temp5 *= 1e-20;
    }
    else if(4.8 < en && en <= 15.0)
    {
      temp1 = en*en;
      temp2 = en*temp1;

      temp5 = 6.113e-4*temp2 -2.11e-2*temp1 +.2216*en -.638556;
      temp5 *= 1e-20;
    }
    
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // dissociation  (E LOSS = 6.0 eV)    
  k=10;
  sigma.at(k).type = COLL_DISS;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).process = "e + O2 -> O(3P) + O(3P) + e  dissociation (6.0 eV)";
  sigma.at(k).E_threshold = 6.0;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5= 0.0;
    if(6.0<= en && en <= 20.)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;

      temp5 = -1.1894e-4*temp3 +6.8655e-3*temp2 -.143425*temp1 +1.26276*en -3.7338513;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(20. < en && en <= 100.)
    {
      temp1 = en*en;

      temp5 = 9.9341e-6*temp1 -1.7857e-3*en +7.924e-2;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // dissociation  (E LOSS = 8.4 eV)    
  k=11;
  sigma.at(k).type = COLL_DISS;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).process = "e + O2 -> O(3P) + O(1D) + e  dissociation (8.4 eV)";
  sigma.at(k).E_threshold = 8.4;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5= 0.0;
    if(8.4<= en && en <= 9.4)
    {
      temp5 = en -8.4;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(9.4 < en && en <= 100.)
    {
      temp1 = en*en;

      temp5 = -1.08852e-4*temp1 +1.10145e-2*en +.92302246;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(100. < en) temp5 = 9.4e-21;
  
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // dissociation  (E LOSS = 9.97 eV)    
  k=12;
  sigma.at(k).type = COLL_DISS;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).process = "e + O2 -> O(1D) + O(1D) + e  dissociation (9.97 eV)";
  sigma.at(k).E_threshold = 9.97;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5= 0.0;
    if(10. <= en && en <= 100.)
    {
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;

      temp5 = 1.1656e-9*temp3 -3.1555e-7*temp2 +1.8544e-5*temp1 +9.464e-4*en -.0110422;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(100. < en) temp5= 7e-22;
    
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // ionization  (E LOSS = 12.06 eV)
  k=13;
  sigma.at(k).type = COLL_ION;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = ELE;
  sigma.at(k).product[1] = O2_P_ION;
  sigma.at(k).process = "e + O2 -> e + e + O2+        ionization";
  sigma.at(k).E_threshold = 12.06;    
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);

    temp5 = 0.0;
    if(12.06 <= en && en <= 100.){
      temp1 = en*en;
      temp2 = en*temp1;
      temp3 = en*temp2;

      temp5 = 8.035e-8*temp3 -1.594e-5*temp2 +5.1392e-4*temp1 +.0658*en -.89892;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5 =0.0;
    }
    else if(100. < en) temp5 = 2.9e-20;

    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;
  }
  // overwrite data
  // read Gudmundsson's e-impact ionization c.s.
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/e_impact_ionization_O2_12.06eV.dat_Gud.dat",cs_sets,false,13);


  // dissociative excitation  (E LOSS = 14.7 eV)   130 nm line excitation  
  k=14;
  sigma.at(k).type = COLL_DISS_EXC;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).process = "e + O2 -> e + O + O(3p3P)    dissociative excitation (14.7 eV)";
  sigma.at(k).E_threshold = 14.7;   
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    temp5= 0.0;
    if(14.7<= en && en <= 100.)
    {
      temp1 = en*en;
      temp2 = en*temp1;

      temp5 = 1.21665e-7*temp2 -3.0483e-5*temp1 +2.51713e-3*en -.030335;
      temp5 *= 1e-20;
      if(temp5 <= 0.0) temp5=0.0;
    }
    else if(100. < en) temp5= 3.8e-22;
    
    cs_value = temp5;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // electron impact detachment 
  k=15;
  sigma.at(k).type = COLL_DET;
  sigma.at(k).target = O_M_ION;
  sigma.at(k).product[0] = ELE;
  sigma.at(k).product[1] = O_GAS;
  sigma.at(k).process = "e + O- -> O + 2e             electron impact detachment";
  sigma.at(k).E_threshold = 1.465;  
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    if(en < 1.465) {detach1[i]= 0.0;}
      else           {detach1[i]= -1.0353e-23*en*en +1.152e-21*en +4.113e-20;}
    if(en <97.1)   {detach2[i]= -1.0353e-23*en*en +1.152e-21*en +4.113e-20;}
        else           {detach2[i]= 1.175e-18*log(en)/en;}
    if(en < 30.)   {sigma1=detach1[i];}
    else	             {sigma1=detach2[i];}

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // overwrite data
  // read Gudmundsson's e-impact detachment c.s.
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/El_im_detachment_Gud.txt",cs_sets,false,15);  // process 16

  // dissociative recombination 
  k=16;
  sigma.at(k).type = COLL_DISS_REC;
  sigma.at(k).target = O2_P_ION;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).process = "e + O2+ -> O(3P) + O(1D)     dissociative recombination";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);
    
    sigma1=1.8475e-19*pow(10.0, -.4524*en);

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // overwrite data
  // read Gudmundsson's dissociative recombination c.s.
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/DR_Gud.txt",cs_sets,false,16);   // process 17

  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }
}


//
// fill data with cross sections for O2+ + O2
//  
void fill_cs_O2_O2p(int& sets){
  
  const double REDUCED_MASS_RATIO = 2.0;
  int    k;
  double en, cs_value, sigma1;
 

  comments           = "CS set for O2p + O2 ";
  output_name        = "cs_O2_O2p.bin";

  sets = 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
  }

  // Charge exchange
  k=0;
  sigma.at(k).type = COLL_BACK;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "O2+ + O2 -> O2 + O2+         charge echange";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);               
    
    if (en <=10.0) {sigma1=1e-18;} 
    else {sigma1=1e-19 +9e-18/en;}

    cs_value= sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // owervrite data
  // read Gudmundsson's charge exchange c.s. 
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  // note1: this is lab energy (target at rest)
  // note2: additional isotropic scattering with 50% of chx --->>> 75% -->> 50% 
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/O2+_charge_exchange_O2.dat_Gud.dat",cs_sets,false,0,REDUCED_MASS_RATIO); // process 18

  // Isotropic elastic scattering
  k=1;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "O2+ + O2 -> O2+ + O2         isotropic elastic scattering";
  sigma.at(k).E_threshold = 0; 
  
  for(int i=0; i<CS_RANGES; i++){
    sigma.at(k).values.at(i) = sigma.at(0).values.at(i) * 0.5;
  }


  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }

}

//
// fill data with cross sections for O- + O2 /O2+ / O2(a1Dg)
//  
void fill_cs_O2_Om(int& sets){   
  
  const double REDUCED_MASS_RATIO = 1.5;
  int    k;
  double en, cs_value, sigma1;
  

  comments           = "CS set for O- + O2 /O2+ / O2(a1Dg) ";
  output_name        = "cs_O2_Om.bin";

  sets = 4;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
  }

  // elastic scattering
  k=0;
  sigma.at(k).type = COLL_ISO ;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "O- + O2 -> O- + O2           elastic scattering";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);

    // some previous formula
    // if (en <=0.2) {sigma1=5e-19;}
    // else	{sigma1=1e-19 +1.79e-19/sqrt(en);}

    // Langevin CS
    sigma1 = 1.6859e-19*sqrt(0.381/en);
  
    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // overwrite data
  // read Gudmundsson's heavy particle cross sections:
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  // note1: this is lab energy (target at rest)
  // interpolation_type="LOG";
  // read_cs_LxCat("cross_sections_O2/Elastic_scat_Om_Gud.txt",cs_sets,false,0,REDUCED_MASS_RATIO); // process 19

  // detachment
  k=1;
  sigma.at(k).type = COLL_DET ;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).product[0] = ELE;
  sigma.at(k).product[1] = O_GAS;
  sigma.at(k).process = "O- + O2 -> O2 + O + e        detachment";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);

    if(en < 1.465) {sigma1 = 0.0;}
    else if((1.465 <= en) && (en < 28.626)) { sigma1 = 2.761e-21*(en - 1.465); }
    else { sigma1 = 7.5e-20; }
  
    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // overwrite data 
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/DirectDet_Gud.txt",cs_sets,false,1); //process 20

  // mutual neutralization
  k=2;
  sigma.at(k).type = COLL_NEUT ;
  sigma.at(k).target = O2_P_ION;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).product[1] = O2_GAS;
  sigma.at(k).process = "O- + O2+ -> O + O2           mutual neutralization";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);

    sigma1=2e-16*pow(10.0, -.0182*en);
  
    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // overwrite data 
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/MneutralizationGud.txt",cs_sets,false,2);// process 21
  

  // associative detachment O- + O2(a1delta)
  k=3;
  sigma.at(k).type = COLL_ASS_DET ;
  sigma.at(k).target = O2_META;
  sigma.at(k).product[0] = ELE;
  sigma.at(k).product[1] = O3_GAS;
  sigma.at(k).process = "O- + O2(a1delta) -> O3 + e   associative detachment";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1) * REDUCED_MASS_RATIO;                // energy conversion: com -> lab;

    // Bronold et al. J. Phys. D 40 6583 (2004) - formula (9)
    // https://doi.org/10.1088/0022-3727/40/21/018

    sigma1 = 5.96/sqrt(en)*1.0e-20; // [m^2]  

    // source: Gudmundsson and Lieberman 2015 Plasma Sources Sci. Technol. 24 035016
    // https://doi.org/10.1088/0963-0252/24/3/035016

    // if (en > 0.134) en = 0.134;
    // sigma1 = sqrt(0.134/en)*5.75e-20; // [m^2]   

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }

}


void fill_cs_Ne_O2p(int& sets){
  
  const double REDUCED_MASS_RATIO = 1.63;    // O2 mass vs. O2+Ne reduced mass
  int    k;
  double en, cs_value, sigma1;
 

  comments           = "CS set for O2p + Ne ";
  output_name        = "cs_Ne_O2p.bin";

  sets = 1;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
  }

  // isotropic elastic scattering with Langevin cross section
  k=0;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).target = NE_GAS;
  sigma.at(k).process = "O2+ + Ne -> O2+ + Ne         isotropic elastic scattering";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);               
    
    //sigma1 = sqrt(0.381e-30*E_CHARGE*PI/2.0/EPSILON0/en[eV]);
    sigma1 = 1.6859e-19*sqrt(0.381/en);

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }

}


//
// fill data with cross sections for O- + Ne / Ne+
//  
void fill_cs_Ne_Om(int& sets){   
  
  const double REDUCED_MASS_RATIO = 1.79;    // O- mass vs. O- + Ne reduced mass
  int    k;
  double en, cs_value, sigma1;
 

  comments           = "CS set for Om + Ne ";
  output_name        = "cs_Ne_Om.bin";

  sets = 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
  }

  // isotropic elastic scattering with Langevin cross section
  k=0;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).target = NE_GAS;
  sigma.at(k).process = "O- + Ne -> O- + Ne         isotropic elastic scattering";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);               
    
    //sigma1 = sqrt(0.381e-30*E_CHARGE*PI/2.0/EPSILON0/en[eV]);
    sigma1 = 1.6859e-19*sqrt(0.381/en);

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // mutual neutralization
  k=1;
  sigma.at(k).type = COLL_NEUT ;
  sigma.at(k).target = NE_ION;
  sigma.at(k).product[0] = O_GAS;
  sigma.at(k).product[1] = NE_GAS;
  sigma.at(k).process = "O- + Ne+ -> O + Ne           mutual neutralization";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);

    sigma1=2e-16*pow(10.0, -.0182*en);
  
    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }
  // overwrite data 
  // source: Gudmundsson et al 2013 Plasma Sources Sci. Technol. 22 035011
  // https://doi.org/10.1088/0963-0252/22/3/035011
  interpolation_type="LOG";
  read_cs_LxCat("cross_sections_O2/MneutralizationGud.txt",cs_sets,false,REDUCED_MASS_RATIO);// process 21


  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }

}

//
// fill data with cross sections for O- + Ne / Ne+
//  
void fill_cs_O2_Nep(int& sets){   
  
  const double REDUCED_MASS_RATIO = 2.58;    // Ne+ mass vs. O2 + Ne reduced mass
  int    k;
  double en, cs_value, sigma1;
 

  comments           = "CS set for Ne+ + O2 ";
  output_name        = "cs_O2_Nep.bin";

  sets = 2;
  cout << "cross section sets found --> "<< sets << endl;
  sigma.resize(sets);
  for (k=0; k<sets; k++){
    sigma.at(k).values.resize(CS_RANGES);
    sigma.at(k).process.reserve(CS_STRING_LENGTH);
    sigma.at(k).product[0] = -1; sigma.at(k).product[1] = -1; sigma.at(k).product[2] = -1; sigma.at(k).product[3] = -1;
  }

  // asymmetric charge transfer process
  // ref: Schlumbohm H 1969 Zeitschrift fur Naturforschung A24 1720–1724
  k=0;
  sigma.at(k).type = COLL_CHARGE_TRANS;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "Ne+ + O2 -> Ne + O2+         charge transfer (one way)";
  sigma.at(k).product[0] = O2_P_ION;
  sigma.at(k).product[1] = NE_GAS;
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);               
    
    
    sigma1 = 3.7e-20;

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }

  // isotropic elastic scattering with Langevin cross section
  k=1;
  sigma.at(k).type = COLL_ISO;
  sigma.at(k).target = O2_GAS;
  sigma.at(k).process = "Ne+ + O2 -> Ne+ + O2         isotropic elastic scattering";
  sigma.at(k).E_threshold = 0; 
  for(int i=0; i<CS_RANGES; i++){
    en = DE_CS * max(i, 1);               
    
    //sigma1 = sqrt(1.562e-30*E_CHARGE*PI/2.0/EPSILON0/en[eV]);
    sigma1 = 1.6859e-19*sqrt(1.562/en);

    cs_value = sigma1;
    sigma.at(k).values.at(i) = (float) cs_value;  
  }


  for (k=0; k<sets; k++){
    cout << sigma.at(k).process << " / type: " << sigma.at(k).type << " / value: " << scientific <<  sigma.at(k).E_threshold << endl;
  }

}


int main() {
  
  interpolation_type.reserve(CS_STRING_LENGTH);     // needed for fixing size in memory
  comments.reserve(CS_STRING_LENGTH);
  input_name.reserve(CS_STRING_LENGTH);
  output_name.reserve(CS_STRING_LENGTH);
  
  printf("PICit! Cross Section Preprocessor\n");
  
  //  process_ele_He_Biagi();
  //  process_He_He_Phelps();
  
  process_ele_Ar_Phelps();
  //  process_ele_Ar_LISBON();
  //  process_ele_Ar_Biagi();
  //  process_ele_Ar_BSR();
  //  process_ele_Ar_COP();
  //  process_ele_Ar_Hayashi();
  //  process_ele_Ar_Morgan();
  //  process_ele_Ar_Puech();
  
  write_cs(output_name, comments);
  
  
  //  process_Arp_Ar_Phelps_website();   // not recommended due to referencing issues
  process_Arp_Ar_Phelps_1994paper();
  write_cs(output_name, comments);
 
  //  process_ele_Ne2p1_excit_Biagi();
  //  write_cs(output_name, comments);
  
  //  process_ele_Ne2p1_excit_Puech();
  //  write_cs(output_name, comments);


  //  process_ele_Ne2p1_excit_BSR();
  //  write_cs(output_name, comments);
  

  //fill_cs_O2_ele(cs_sets);
  //write_cs(output_name, comments);  
 //
  //process_ele_Ne_Biagi71();
  //write_cs(output_name, comments);
//
  //merge_CS_files("cs_NeO2_e.bin", "cs_O2_e.bin", "cs_Ne_e.bin");
  //remove("cs_Ne_e.bin");
  //remove("cs_O2_e.bin");
//
//
  //fill_cs_O2_O2p(cs_sets);
  //write_cs(output_name, comments);  
 //
  //fill_cs_Ne_O2p(cs_sets);
  //write_cs(output_name, comments);  
//
  //merge_CS_files("cs_NeO2_O2p.bin", "cs_O2_O2p.bin", "cs_Ne_O2p.bin");
  //remove("cs_Ne_O2p.bin");
  //remove("cs_O2_O2p.bin");
//
//
//
  //fill_cs_O2_Om(cs_sets);
  //write_cs(output_name, comments);  
//
  //fill_cs_Ne_Om(cs_sets);
  //write_cs(output_name, comments);  
//
  //merge_CS_files("cs_NeO2_Om.bin", "cs_O2_Om.bin", "cs_Ne_Om.bin");
  //remove("cs_Ne_Om.bin");
  //remove("cs_O2_Om.bin");
//
//
//
  //process_Nep_Ne_Jovanovic();
  //write_cs(output_name, comments);
//
  //fill_cs_O2_Nep(cs_sets);
  //write_cs(output_name, comments);
//
  //merge_CS_files("cs_NeO2_Nep.bin", "cs_O2_Nep.bin", "cs_Ne_Nep.bin");
  //remove("cs_Ne_Nep.bin");
  //remove("cs_O2_Nep.bin");

  process_Arp_Ar_Phelps_1994paper();
  write_cs(output_name, comments);

  return 0;
}

