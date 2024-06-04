/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Dallas R. Trinkle (UIUC), dtrinkle@illinois.edu
------------------------------------------------------------------------- */

#include "fix_vsgcmc.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixVirtualSemiGrandCanonicalMC::FixVirtualSemiGrandCanonicalMC(class LAMMPS_NS::LAMMPS *, int narg, char **arg) :
        Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), type_list(nullptr),
        qtype(nullptr), local_swap_atom_list(nullptr), random_equal(nullptr),
        random_unequal(nullptr), c_pe(nullptr)
{
  if (narg < 10) error->all(FLERR, "Illegal fix vsgcmc command");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  ncycles = utils::inumeric(FLERR, arg[4], false, lmp);
  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  double temperature = utils::numeric(FLERR, arg[6], false, lmp);

  if (nevery <= 0) error->all(FLERR, "Illegal fix vsgcmc command");
  if (ncycles < 0) error->all(FLERR, "Illegal fix vsgcmc command");
  if (seed <= 0) error->all(FLERR, "Illegal fix vsgcmc command");
  if (temperature <= 0.0) error->all(FLERR, "Illegal fix vsgcmc command");

  beta = 1.0 / (force->boltz * temperature);

  memory->create(type_list, atom->ntypes, "vsgcmc:type_list");

  // read options from end of input line

  options(narg - 7, &arg[7]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp, seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp, seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  mc_active = 0;


  atom_swap_nmax = 0;
  local_swap_atom_list = nullptr;

  // set comm size needed by this Fix

  if (atom->q_flag)
    comm_forward = 2;
  else
    comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixVirtualSemiGrandCanonicalMC::~FixVirtualSemiGrandCanonicalMC()
{
  memory->destroy(type_list);
  memory->destroy(list_type);
  memory->destroy(qtype);
  memory->destroy(chemdifferences);
  memory->destroy(swapchem);
  memory->destroy(swapindex);
  memory->destroy(nattempt);
  memory->destroy(chempotave);
  delete random_equal;
  delete random_unequal;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal fix vsgcmc command");

  nswaptypes = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "types") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix vsgcmc command");
      iarg++;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        if (nswaptypes >= atom->ntypes) error->all(FLERR, "Illegal fix vsgcmc command");
        type_list[nswaptypes] = utils::numeric(FLERR, arg[iarg], false, lmp);
        nswaptypes++;
        iarg++;
      }
    } else
      error->all(FLERR, "Illegal fix vsgcmc command");
  }
}

/* ---------------------------------------------------------------------- */

int FixVirtualSemiGrandCanonicalMC::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::init()
{
  if (!atom->mass) error->all(FLERR, "Fix vsgcmc requires per atom type masses");
  if (atom->rmass_flag && (comm->me == 0))
    error->warning(FLERR, "Fix vsgcmc will use per-type masses for velocity rescaling");

  c_pe = modify->get_compute_by_id("thermo_pe");

  int *type = atom->type;

  if (nswaptypes < 2) error->all(FLERR, "Must specify at least 2 types in fix vsgcmc command");

  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    if (type_list[iswaptype] <= 0 || type_list[iswaptype] > atom->ntypes)
      error->all(FLERR, "Invalid atom type in fix vsgcmc command");

  memory->create(list_type, atom->ntypes+1, "vsgmcm:list_type");
  for (int i=0; i<=atom->ntypes; ++i)
      list_type[i] = 0;
  for (int i=0; i<nswaptypes; ++i) {
      list_type[type_list[i]] = i;
  }

  // now we need to set up some helper arrays to keep track of what we swap, etc.
  // our chemical potential differences are:
  // type[1]-type[0], type[2]-type[0], ..., type[nswap-1]-type[0], type[2]-type[1], ...
  // ... type[nswap-1]-type[nswap-2]
  // For type[1]-type[0], we need exp(-beta*dE(0->1)) and exp(+beta*dE(1->0)),
  // but we have to keep those two averages separate.
  // swapchem[i]: list of indices to calculate swaps = (0,1,...,i-1,i+1,...,N-1)
  // swapindex[i]: list of which chemical potential difference that is
  // chemdifferences[d]: [0] is the + and [1] is the - (A-B+ == A->B)
  nchempot = (nswaptypes*(nswaptypes-1));
  memory->create(swapchem, nswaptypes, nswaptypes-1, "vsgcmc:swapchem");
  memory->create(swapindex, nswaptypes, nswaptypes-1, "vsgcmc:swapindex");
  memory->create(chemdifferences, nchempot, 2, "vsgcmc:chemdifferences");

  int nchem = 0;
  for (int i=0; i<nswaptypes; ++i) {
      for (int j=i+1; j<nswaptypes; ++j) {
          chemdifferences[nchem][0] = j;
          chemdifferences[nchem][1] = i;
          swapchem[i][j-1] = j;
          swapindex[i][j-1] = nchem;
          nchem++;
          chemdifferences[nchem][0] = i;
          chemdifferences[nchem][1] = j;
          swapchem[j][i] = i;
          swapindex[j][i] = nchem;
          nchem++;
      }
  }

  memory->create(nattempt, nchempot, "vsgcmc:nattempt");
  memory->create(chempotave, nchempot, "vsgcmc:chempotave");
  for (int i=0; i<nchempot; ++i) {
      nattempt[i] = 0;
      chempotave[i] = 0.;
  }

  // this is only required for non-semi-grand
  // in which case, nswaptypes = 2

  if (atom->q_flag) {
    double qmax, qmin;
    int firstall, first;
    memory->create(qtype, nswaptypes, "vsgcmc:qtype");
    for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++) {
      first = 1;
      for (int i = 0; i < atom->nlocal; i++) {
        if (atom->mask[i] & groupbit) {
          if (type[i] == type_list[iswaptype]) {
            if (first) {
              qtype[iswaptype] = atom->q[i];
              first = 0;
            } else if (qtype[iswaptype] != atom->q[i])
              error->one(FLERR, "All atoms of a swapped type must have the same charge.");
          }
        }
      }
      MPI_Allreduce(&first, &firstall, 1, MPI_INT, MPI_MIN, world);
      if (firstall)
        error->all(FLERR,
                   "At least one atom of each swapped type must be present to define charges.");
      if (first) qtype[iswaptype] = -DBL_MAX;
      MPI_Allreduce(&qtype[iswaptype], &qmax, 1, MPI_DOUBLE, MPI_MAX, world);
      if (first) qtype[iswaptype] = DBL_MAX;
      MPI_Allreduce(&qtype[iswaptype], &qmin, 1, MPI_DOUBLE, MPI_MIN, world);
      if (qmax != qmin) error->all(FLERR, "All atoms of a swapped type must have same charge.");
    }
  }

  // check to see if itype and jtype cutoffs are the same
  // if not, reneighboring will be needed between swaps

  double **cutsq = force->pair->cutsq;
  unequal_cutoffs = false;
  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    for (int jswaptype = 0; jswaptype < nswaptypes; jswaptype++)
      for (int ktype = 1; ktype <= atom->ntypes; ktype++)
        if (cutsq[type_list[iswaptype]][ktype] != cutsq[type_list[jswaptype]][ktype])
          unequal_cutoffs = true;

  // check that no swappable atoms are in atom->firstgroup
  // swapping such an atom might not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);

    if (flagall) error->all(FLERR, "Cannot do vsgcmc on atoms in atom_modify first group");
  }
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo swaps
------------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  mc_active = 1;

  // ensure current system is ready to compute energy

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);

  // energy_stored = energy of current state
  // will be updated after accepted swaps

  energy_stored = energy_full();

  // attempt Ncycle atom swaps

  update_atoms_list();
  for (int i = 0; i < ncycles; i++) virtual_semi_grand();

  // update time counter
  next_reneighbor = update->ntimestep + nevery;

  mc_active = 0;
}

/* ----------------------------------------------------------------------
   calculate a virtual semd-grand exchange of a single atom, store
   exp ( -beta*dE )
   in the appropriate average
   NOTE: atom charges are assumed equal and so are not updated
------------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::virtual_semi_grand()
{
  if (nswap == 0) return 0;

  // pre-swap energy

  double energy_before = energy_stored;

  // pick a random atom and perform all transmutations on it

  int itype, jtype, jswaptype, nchem, i_ind;
  int i = pick_semi_grand_atom();
  if (i >= 0) {
      itype = atom->type[i];
      i_ind = list_type[itype];
  }
//  if (i >= 0) {
//    jswaptype = static_cast<int>(nswaptypes * random_unequal->uniform());
//    jtype = type_list[jswaptype];
//    itype = atom->type[i];
//    while (itype == jtype) {
//      jswaptype = static_cast<int>(nswaptypes * random_unequal->uniform());
//      jtype = type_list[jswaptype];
//    }
//    atom->type[i] = jtype;
//  }

  // if unequal_cutoffs, call comm->borders() and rebuild neighbor list
  // else communicate ghost atoms
  // call to comm->exchange() is a no-op but clears ghost atoms
  for (int j_ind=0; j_ind<nswaptypes-1; ++j_ind) {
      if (i >= 0) {
          jswaptype = type_list[swapchem[i_ind][j_ind]];
          nchem = swapindex[i_ind][j_ind];
          atom->type[i] = jswaptype;
      }
      if (unequal_cutoffs) {
          if (domain->triclinic) domain->x2lamda(atom->nlocal);
          comm->exchange();
          comm->borders();
          if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
          if (modify->n_pre_neighbor) modify->pre_neighbor();
          neighbor->build(1);
      } else {
          comm->forward_comm(this);
      }

      // post-swap energy

      if (force->kspace) force->kspace->qsum_qsq();
      double energy_after = energy_full();

      // now to store in the appropriate average:
      nattempt[nchem]++;
      double incr_chem_pot = exp(beta * (energy_before - energy_after)) - chempotave[nchem];
      chempotave[nchem] += incr_chem_pot/nattempt[nchem];
  }

  // restore the swapped atom
  // do not need to re-call comm->borders() and rebuild neighbor list
  //   since will be done on next cycle or in Verlet when this fix finishes

  if (i >= 0) atom->type[i] = itype;
  if (force->kspace) force->kspace->qsum_qsq();

  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixVirtualSemiGrandCanonicalMC::energy_full()
{
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  if (modify->n_post_force_any) modify->post_force(vflag);

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixVirtualSemiGrandCanonicalMC::pick_semi_grand_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int>(nswap * random_equal->uniform());
  if ((iwhichglobal >= nswap_before) && (iwhichglobal < nswap_before + nswap_local)) {
    int iwhichlocal = iwhichglobal - nswap_before;
    i = local_swap_atom_list[iwhichlocal];
  }

  return i;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::update_atoms_list()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  if (atom->nmax > atom_swap_nmax) {
    memory->sfree(local_swap_atom_list);
    atom_swap_nmax = atom->nmax;
    local_swap_atom_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_atom_list");
  }

  nswap_local = 0;

  if (region) {
    for (int i = 0; i < nlocal; i++) {
      if (region->match(x[i][0], x[i][1], x[i][2]) == 1) {
        if (atom->mask[i] & groupbit) {
          int itype = atom->type[i];
          int iswaptype;
          for (iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
            if (itype == type_list[iswaptype]) break;
          if (iswaptype == nswaptypes) continue;
          local_swap_atom_list[nswap_local] = i;
          nswap_local++;
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        int itype = atom->type[i];
        int iswaptype;
        for (iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
          if (itype == type_list[iswaptype]) break;
        if (iswaptype == nswaptypes) continue;
        local_swap_atom_list[nswap_local] = i;
        nswap_local++;
      }
    }
  }

  MPI_Allreduce(&nswap_local, &nswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&nswap_local, &nswap_before, 1, MPI_INT, MPI_SUM, world);
  nswap_before -= nswap_local;
}

/* ---------------------------------------------------------------------- */

int FixVirtualSemiGrandCanonicalMC::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;

  if (atom->q_flag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
      buf[m++] = q[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;
  last = first + n;

  if (atom->q_flag) {
    for (i = first; i < last; i++) {
      type[i] = static_cast<int>(buf[m++]);
      q[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++) type[i] = static_cast<int>(buf[m++]);
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratio
------------------------------------------------------------------------- */

double FixVirtualSemiGrandCanonicalMC::compute_vector(int n)
{
  if (n == 0) return nswap_attempts;
  if (n == 1) return nswap_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixVirtualSemiGrandCanonicalMC::memory_usage()
{
  double bytes = (double) atom_swap_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = nswap_attempts;
  list[n++] = nswap_successes;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), n, fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixVirtualSemiGrandCanonicalMC::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;

  seed = static_cast<int>(list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int>(list[n++]);
  random_unequal->reset(seed);

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR, "Must not reset timestep when restarting fix vsgcmc");
}

/* ----------------------------------------------------------------------
   extract variable which stores whether MC is active or not
     active = MC moves are taking place
     not active = normal MD is taking place
------------------------------------------------------------------------- */

void *FixVirtualSemiGrandCanonicalMC::extract(const char *name, int &dim)
{
  if (strcmp(name,"mc_active") == 0) {
    dim = 0;
    return (void *) &mc_active;
  }
  return nullptr;
}
