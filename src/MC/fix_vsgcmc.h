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
  void *extract(const char *, int &) override;

 private:
  int nevery, seed;
  int ncycles;
  int niswap, njswap;                  // # of i,j swap atoms on all procs
  int niswap_local, njswap_local;      // # of swap atoms on this proc
  int niswap_before, njswap_before;    // # of swap atoms on procs < this proc
  int nswap;                           // # of swap atoms on all procs
  int nswap_local;                     // # of swap atoms on this proc
  int nswap_before;                    // # of swap atoms on procs < this proc
  class Region *region;                // swap region
  char *idregion;                      // swap region id

  int mc_active;              // 1 during MC trials, otherwise 0

  int nswaptypes;
  int *type_list;

  int nchempot;
  int **swapindex;
  int **swapchem;
  int **swapsign;
  int **chemdifferences;

  double *mu;

  double nswap_attempts;
  double nswap_successes;

  bool unequal_cutoffs;

  int atom_swap_nmax;
  double beta;
  double *qtype;
  double energy_stored;
  int *local_swap_iatom_list;
  int *local_swap_jatom_list;
  int *local_swap_atom_list;

  class RanPark *random_equal;
  class RanPark *random_unequal;

  class Compute *c_pe;

  void options(int, char **);
  int attempt_virtual_semi_grand();
  double energy_full();
  int pick_semi_grand_atom();
  int pick_i_swap_atom();
  int pick_j_swap_atom();
};

}    // namespace LAMMPS_NS

#endif
#endif
