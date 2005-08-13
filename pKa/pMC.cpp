//
// $Id$
//
// Jens Erik Nielsen, University College Dublin 2005
//
#include "pMC.h"
//
//
//

void MC::reformat_arrays() {
  //
  // Make the matrix
  //
  _groups=static_cast<int>(_intpKas.size());
  int count=0;
  for (int row=0;row<_groups;row++) {
    vector<double> row_vals;
    for (int column=0;column<_groups;column++) {
      row_vals.push_back(_lin_matrix[count]);
      count=count+1;
    }
    _matrix.push_back(row_vals);
  }
  //
  // Set natural log
  //
  lnten=log(10);
  return;
}

//
// ---------------------
//

vector<double> MC::calc_pKa(double pH_start,double pH_end, double pH_step) {
  //
  // Calculate pKa values for the system
  //
  // First get charges at all pH values
  //
  vector< vector<double> > charges;
  for (double pH=pH_start;pH<pH_end;pH=pH+pH_step) {
    //printf ("%5.3f\n",pH);
    charges.push_back(calc_charge(pH));
  }
  // 
  // Now determine pKa values
  //
  vector<double> pKas;
  for (int group=0;group<_groups;group++) {
     int count=0;
    double pKa=-9999.9;
    double last_crg=charges[count][group];
    for (double pH=pH_start;pH<pH_end;pH=pH+pH_step) {
      double this_crg=charges[count][group];
      if (_acid_base[group]==1) {
	if (this_crg<=0.5 && last_crg>0.5) {
	  pKa=(last_crg-0.5)/(last_crg-this_crg)*pH_step+(pH-pH_step);
	}
      } else {
	if (this_crg<=-0.5 && last_crg>-0.5) {
	  pKa=(last_crg-(-0.5))/(last_crg-this_crg)*pH_step+(pH-pH_step);
	}
      }
      last_crg=this_crg;
      count=count+1;
    }
    pKas.push_back(pKa);
  }   
  return pKas;
}

//
// ---------------------
//

vector<double> MC::calc_charge(double pH) {
  //
  // Calculate the fractional charges at this pH
  //
  // Initialise random number generator
  //
  srand(time(NULL));
  //
  // Get a random starting state
  //
  vector<int> current_state;
  vector<int> try_state;
  vector<int> sum_state;
  for (int group=0;group<_groups;group++) {
    current_state.push_back(static_cast<int>(rand()%2));
    if (current_state[group]==2) {
      current_state[group]=1;
    }
    // 
    // Dummy initialisation of try_state and sum_state
    //
    try_state.push_back(0);
    sum_state.push_back(0);
  }
  //
  // 
  //
  float current_energy=get_energy(pH,current_state);
  //
  // Start the MC loop
  //
  int eqsteps=_MCsteps/10;
  int keep=0;
  double tilf=0.0;
  for (int step=0;step<_MCsteps;step++) {
    //
    // Copy the current state to trystate
    //
    for (int count=0;count<_groups;count++) {
      try_state[count]=current_state[count];
    }
    //
    // Change a random group
    //
    int rand_group=rand()%_groups;
    try_state[rand_group]=abs(try_state[rand_group]-1);
    //
    // Get the energy of the new state
    //
    float try_energy=get_energy(pH,try_state);
    //
    // Keep or reject?
    //
    float diff=try_energy-current_energy;
    keep=0;
    if (diff<0.0) {
      keep=1;
    } else {
      if (diff<20.0) {
	tilf=static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)+static_cast<double>(1));
	if (tilf<exp(-diff)) keep=1;
      }
    }
    //
    // If we keep the state then copy to current state
    //
    if (keep==1) {
      for (int count=0;count<_groups;count++) {
	current_state[count]=try_state[count];
      }
      current_energy=try_energy;
    }
    //
    // Record the state if we have equilibrated
    //
    if (step>eqsteps) {
      for (int count=0;count<_groups;count++) {
	sum_state[count]=sum_state[count]+current_state[count];
      }
    }
  }
  //
  // Calculate fractional charge
  //
  vector<double> charges_thispH;
  for (int count=0;count<_groups;count++) {
    charges_thispH.push_back(static_cast<double>(_acid_base[count])\
			     *static_cast<double>(sum_state[count])/\
			     static_cast<double>(_MCsteps-eqsteps));
  }
  return charges_thispH;
}

//
// --------------------
// 

double MC::get_energy(double pH,vector<int> state) {
  //
  // Calculate the energy of the present state
  //
  double energy=0.0;
  for (int group1=0;group1<_groups;group1++) {
    //
    // Add the energy from the intrinsic pKa
    //
    if (state[group1]==1) {
      energy=energy+_acid_base[group1]*lnten*(pH-_intpKas[group1]);
      //
      // Add the charged-charged energies
      //
      for (int group2=0;group2<_groups;group2++) {
	if (state[group2]==1 && group2!=group1) {
	  energy=energy+_matrix[group1][group2]/2.0;
	}
      }
    }
  }
  return energy;
}
	
  
