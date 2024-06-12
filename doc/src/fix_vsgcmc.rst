.. index:: fix vsgcmc

fix vsgcmc command
==================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID vsgcmc N T keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* vsgcmc = style name of this fix command
* N = invoke this fix every N steps
* T = temperature of the system (temperature units)
* keyword = *types* or *region*

  .. parsed-literal::

     keyword = *types* or *region*
       *types* values = two or more atom types for transmutations (required)
       *region* value = region-ID
         region-ID = ID of region to use as a transmutation volume

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all vsgcmc 100 300 types 1 2
   fix 2 surface vsgcmc 20 1000 types 2 3 4
   fix 3 all vsgcmc 100 1200 region interactionlayer types 2 4

Description
"""""""""""

This fix performs semi-grand canonical Widom method (or virtual semi-grand
canonical Monte Carlo) exchange of atomic types at the given
temperature as derived in :ref:`(Sindzingre) <Sindzingre1>`, and further
described in :ref:`(Anwar) <Anwar1>`. This allows the computation of the excess
chemical potential difference between two (or more) chemical species.

Every N timesteps the fix sweeps through all atoms in the group (or region);
for each atom of an identified type, it transmutes it into each of the other
types, recording the energy difference. The Boltzmann factor of energy
change is averaged, and reported.

The *types* keyword is required. At least two atom types must be
specified.

This command may optionally use the *region* keyword to define transmuting
volume.  The specified region must have been previously defined with a
:doc:`region <region>` command.  It must be defined with side = *in*\ .
Only atoms in the specified region are transmuted.

Note that neighbor lists are re-built every timestep that this fix is
invoked, so you should not set N to be too small. As the computation
of energy change involves recomputing the total energy of the entire
cell, it is recommended to not let N be smaller than the number of atoms.
See the :doc:`neighbor <neighbor>` command for details.

The difference in excess chemical potential dmu_ex is defined as:

.. math::

   \delta\mu_{ex} = \mu_{ex,B} - \mu_{ex,A} = -kT \ln(<\exp(-\Delta U(A-B+)/{k_B T}>)

where :math:`k_B` is the Boltzmann constant, :math:`T` is the
user-specified temperature, :math:`\Delta U(A-B+)` is the potential energy
change from transmuting one atom of type A into type B.

This fix computes a global vector of length Ntypes*(Ntypes-1), which can
be accessed by various :doc:`output commands <Howto_output>`.  The vector
values are the following quantities, recomputed every N steps:

* 1 = :math:`<\exp(-\Delta U(t_1\to t_2)/{k T}>)`
* 2 = :math:`<\exp(-\Delta U(t_2\to t_1)/{k T}>)`
* 3 = :math:`<\exp(-\Delta U(t_1\to t_3)/{k T}>)`
* 4 = :math:`<\exp(-\Delta U(t_3\to t_1)/{k T}>)`
* ...
* Ntypes*(Ntypes-1)-1 = :math:`<\exp(-\Delta U(t_{Ntypes-1}\to t_{Ntypes})/{k T}>)`
* Ntypes*(Ntypes-1) = :math:`<\exp(-\Delta U(t_{Ntypes}\to t_{Ntypes-1})/{k T}>)`

following the order of the types listed with the *types* keyword
The vector values calculated by this fix are "intensive". When this fix is initialized,
the indices and their corresponding transmutation pairs are written to the logfile.
Each time the fix is called, the averages are computed over all atoms of each type.

Some fixes have an associated potential energy. Examples of such fixes
include: :doc:`efield <fix_efield>`, :doc:`gravity <fix_gravity>`,
:doc:`addforce <fix_addforce>`, :doc:`restrain <fix_restrain>`, and
:doc:`wall fixes <fix_wall>`.  For that energy to be included in the
total potential energy of the system (the quantity used when performing
exchanges), you MUST enable the :doc:`fix_modify <fix_modify>`
*energy* option for that fix.  The doc pages for individual :doc:`fix
<fix>` commands specify if this should be done.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the fix to :doc:`binary restart files
<restart>`.  This includes information about the next timestep for
transmutations. See the :doc:`read_restart <read_restart>` command for
info on how to re-specify a fix in an input script that reads a restart
file, so that the operation of the fix continues in an uninterrupted fashion.

.. note::

   For this to work correctly, the timestep must **not** be changed
   after reading the restart with :doc:`reset_timestep
   <reset_timestep>`.  The fix will try to detect it and stop with an
   error.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Do not set "neigh_modify once yes" or else this fix will never be
called.  Reneighboring is **required**.

This fix style requires an :doc:`atom style <atom_style>` with per atom
type masses.

Can be run in parallel, but some aspects of the transmutation procedure
will not scale well in parallel. Only usable for 3D simulations.


Related commands
""""""""""""""""

:doc:`fix atom/swap <fix_atom_swap>`,
:doc:`fix widom <fix_widom>`,
:doc:`neighbor <neighbor>`


Default
"""""""

There are no defaults.

----------

.. _Sindzingre1:

**(Sindzingre)** P. Sindzingre, G. Ciccotti, C. Massobrio, and D. Frenkel,
"Partial enthalpies and related quantities in mixtures from computer simulation."
*Chem. Phys. Lett.* **136**, 35-41 (1987). doi:10.1016/0009-2614(87)87294-9

.. _Anwar1:

**(Anwar)** J. Anwar, C. Leitold and B. Peters,
"Solid-solid phase equilibria in the NaCl-KCl system."
*J. Chem. Phys.* **152**, 144109 (2020). doi: 10.1063/5.0003224
