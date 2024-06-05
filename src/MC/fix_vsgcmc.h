/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------

 Implementation of the semigrand canonical Widom method, from
 P. Sindzingre, G. Ciccotti, C. Massobrio, and D. Frenkel,
 Chem. Phys. Lett. 136, 35 (1987); applied in
 J. Anwar, C. Leitold and B. Peters, J. Chem. Phys. 152, 144109 (2020)
 doi: 10.1063/5.0003224

 It calculates the average of exp(-beta*dE(A-,B+)) during a simulation;
 this quantity can be then related to the difference in the excess
 chemical potential between species B and A.

 2024 June: Dallas R. Trinkle
 */

#ifdef FIX_CLASS
// clang-format off
FixStyle(vsgcmc,FixVirtualSemiGrandCanonicalMC);
// clang-format on
#else

#ifndef LMP_FIX_VIRTUAL_SEMIGRANDCANONICAL_H
#define LMP_FIX_VIRTUAL_SEMIGRANDCANONICAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixVirtualSemiGrandCanonicalMC : public Fix {
 public:
  FixVirtualSemiGrandCanonicalMC(class LAMMPS *, int, char **);
  ~FixVirtualSemiGrandCanonicalMC() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double compute_vector(int) override;
  double memory_usage() override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  int nevery;
//  int seed;
//  int ncycles;
  int nswap;                           // # of swap atoms on all procs
  int nswap_local;                     // # of swap atoms on this proc
  int nswap_before;                    // # of swap atoms on procs < this proc
  class Region *region;                // swap region
  char *idregion;                      // swap region id

  int nswaptypes;
  int *type_list;
  int *list_type; // inverse mapping of type_list

  int nchempot;
  int **swapindex;
  int **swapchem;
  int **chemdifferences;

  int *nattempt;
  double *chempotave;

  bool unequal_cutoffs;

  int atom_swap_nmax;
  double beta;
  double *qtype;
  double energy_stored;
  int *local_swap_atom_list;

  class Compute *c_pe;

  void options(int, char **);
  void virtual_semi_grand(int);
  double energy_full();
  int pick_semi_grand_atom(int);
  void update_atoms_list();
};

}    // namespace LAMMPS_NS

#endif
#endif
