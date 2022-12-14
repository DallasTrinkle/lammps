/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NSTENCIL_CLASS
// clang-format off
NStencilStyle(full/multi/2d,
              NStencilFullMulti2d, NS_FULL | NS_MULTI | NS_2D | NS_ORTHO | NS_TRI);
// clang-format on
#else

#ifndef LMP_NSTENCIL_FULL_MULTI_2D_H
#define LMP_NSTENCIL_FULL_MULTI_2D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilFullMulti2d : public NStencil {
 public:
  NStencilFullMulti2d(class LAMMPS *);
  void create() override;

 protected:
  void set_stencil_properties() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
