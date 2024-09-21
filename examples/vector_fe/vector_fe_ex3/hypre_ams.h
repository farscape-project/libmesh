// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef HYPRE_AMS_H
#define HYPRE_AMS_H

#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fem_system.h"
#include "libmesh/newton_solver.h"
#include "libmesh/petsc_linear_solver.h"

using namespace libMesh;

class HypreAMS
{
  private:
    const FEMSystem & nedelec_sys;
    System & lagrange_sys;
    const MeshBase & mesh;
    const PC & pc;

    dof_id_type n_glb_edges, n_loc_edges, n_glb_vertices, n_loc_vertices;

    PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

    Mat G;
    Vec x, y, z;
    Vec Gx, Gy, Gz;

    void Size();
    void Create();
    void Build();
    void Mult();
    void Set();

  public:
    HypreAMS(FEMSystem & sys);
};

#endif // HYPRE_AMS_H