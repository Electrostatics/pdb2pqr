//
// $Id$
//
// Jens Erik Nielsen, University College Dublin 2005
//
#ifndef PKAMC_H
#define PKAMC_H

#include <string>
#include <vector>
#include <math.h>
#include <time.h>

using namespace std;

class MC {
public:
  //
  // Functions
  //
  MC(vector<double> intpKas, vector<double> lin_matrix,vector<double> ab): _intpKas(intpKas),_lin_matrix(lin_matrix),_acid_base(ab) {
    reformat_arrays();
    // Set default value for MCsteps
    _MCsteps=20000;
  };
  //
  //void set_acid_base(vector<int> acid_base);
  //
  // Reformats the matrix
  void reformat_arrays();
  //
  void set_MCsteps(int MCsteps) {
    _MCsteps=MCsteps;
    return;
  }
  vector<double> calc_pKa(double pH_start,double pH_end, double pH_step);
  vector<double> calc_charge(double);
  double get_energy(double pH, vector<int> state);
  //
  //
  // Variables
  //
  vector<double> _intpKas, _lin_matrix, _acid_base;
  vector<vector<double> > _matrix, _charges;
  int _groups, _MCsteps;
  double lnten;
  //vector<vector <double>> charges;
};
#endif
