/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_full_multi2.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairFullMulti2::NPairFullMulti2(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   multi2-type stencil is itype-jtype dependent   
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void NPairFullMulti2::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,jbin,which,ns,imol,iatom,moltemplate;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  int js;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    ibin = atom2bin_multi2[itype][i];
    
    // loop through stencils for all types
    for (jtype = 1; jtype <= atom->ntypes; jtype++) {

      // if same type use own bin
      if(itype == jtype) jbin = ibin;
	  else jbin = coord2bin(x[i], jtype);

      // loop over all atoms in surrounding bins in stencil including self
      // skip i = j
      // use full stencil for all type combinations

      s = stencil_multi2[itype][jtype];
      ns = nstencil_multi2[itype][jtype];
      
      for (k = 0; k < ns; k++) {
	    js = binhead_multi2[jtype][jbin + s[k]];
	    for (j = js; j >= 0; j = bins_multi2[jtype][j]) {
	      if (i == j) continue;
	  
	      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

	      delx = xtmp - x[j][0];
	      dely = ytmp - x[j][1];
	      delz = ztmp - x[j][2];
	      rsq = delx*delx + dely*dely + delz*delz;
          
	      if (rsq <= cutneighsq[itype][jtype]) {
	        if (molecular != Atom::ATOMIC) {
	          if (!moltemplate)
		        which = find_special(special[i],nspecial[i],tag[j]);
	          else if (imol >= 0)
		        which = find_special(onemols[imol]->special[iatom],
		    		     onemols[imol]->nspecial[iatom],
		    		     tag[j]-tagprev);
	          else which = 0;
	          if (which == 0) neighptr[n++] = j;
	          else if (domain->minimum_image_check(delx,dely,delz))
		        neighptr[n++] = j;
	          else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
	        } else neighptr[n++] = j;
	      }
	    }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
  list->gnum = 0;
}
