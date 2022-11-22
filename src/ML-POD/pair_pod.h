/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mlpod,PairPOD);
// clang-format on
#else

#ifndef LMP_PAIR_POD_H
#define LMP_PAIR_POD_H

#include "pair.h"

namespace LAMMPS_NS {

  class PairPOD : public Pair {
    public:
      PairPOD(class LAMMPS *);
      ~PairPOD() override;
      void compute(int, int) override;

      void settings(int, char **) override;
      void coeff(int, char **) override;
      void init_style() override;
      double init_one(int, int) override;
      double memory_usage() override;

      int dim = 3;

      double *gd;         // global linear descriptors
      double *podcoeff;   // POD coefficients
      double *newpodcoeff;// normalized POD coefficients
      double *energycoeff; // energy coefficients
      double *forcecoeff; // force coefficients

      void estimate_tempmemory();
      void free_tempmemory();
      void allocate_tempmemory();

      void free_memory();
      void allocate_memory();

      void lammpsNeighPairs(double **x, int **firstneigh, int *atomtype, int *map, int *numneigh, int i);

  protected:
    int nablockmax=0;     // maximum number of atoms per computation block
    int nij=0;    //  number of atom pairs
    int nijmax=0;  // maximum number of atom pairs
    int szd=0;    // size of tmpmem

    class CPOD *podptr;

    // temporary arrays for computation blocks

    double *tmpmem=NULL;  // temporary memory
    int *typeai=NULL;     // types of atoms I only
    int *numneighsum=NULL;// cumulative sum for an array of numbers of neighbors
    double *rij=NULL;     // (xj - xi) for all pairs (I, J)
    int *idxi=NULL;       // storing linear indices for all pairs (I, J)
    int *ai=NULL;         // IDs of atoms I for all pairs (I, J)
    int *aj=NULL;         // IDs of atoms J for all pairs (I, J)
    int *ti=NULL;         // types of atoms I for all pairs (I, J)
    int *tj=NULL;         // types of atoms J  for all pairs (I, J)

    double **scale;         // for thermodynamic integration
  };

}    // namespace LAMMPS_NS

#endif
#endif
