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

// The data structures, systems and solvers we may use
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fem_system.h"
#include "libmesh/newton_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

using namespace libMesh;

void BuildHypreAMS(FEMSystem & sys)
{
  // The Nédélec system (and associated mesh) we build the preconditioner for
  const FEMSystem & nedelec_sys = sys;
  const MeshBase & mesh = nedelec_sys.get_mesh();

  // Dummy Lagrange system defined over the same mesh so we can enumerate the vertices
  System & lagrange_sys = sys.get_equation_systems().add_system<System>("dummy");
  lagrange_sys.add_variable("dummy");
  lagrange_sys.get_equation_systems().reinit_mesh();

  // Global (i.e. total) and local (i.e. to this processor) number of edges and vertices
  const dof_id_type n_glb_edges = nedelec_sys.n_dofs();
  const dof_id_type n_loc_edges = nedelec_sys.n_local_dofs();
  const dof_id_type n_glb_vertices = lagrange_sys.n_dofs();
  const dof_id_type n_loc_vertices = lagrange_sys.n_local_dofs();

  // Create the discrete grandient matrix, representing the edges in terms of its vertices
  PetscMatrix<Real> G(mesh.comm());
  G.init(n_glb_edges, n_glb_vertices, n_loc_edges, n_loc_vertices, 2, 2);

  // Create vectors for the coordinates of the vertices
  PetscVector<Real> x(mesh.comm(), n_glb_vertices, n_loc_vertices);
  PetscVector<Real> y(mesh.comm(), n_glb_vertices, n_loc_vertices);
  PetscVector<Real> z(mesh.comm(), n_glb_vertices, n_loc_vertices);

  // Create vectors for the mat-vec products, representing constant vector fields in the Nédélec basis
  PetscVector<Real> Gx(mesh.comm(), n_glb_edges, n_loc_edges);
  PetscVector<Real> Gy(mesh.comm(), n_glb_edges, n_loc_edges);
  PetscVector<Real> Gz(mesh.comm(), n_glb_edges, n_loc_edges);

  // Populate the discrete gradient matrix and the coordinate vectors
  for (const auto & elem : mesh.active_local_element_ptr_range())
    for(unsigned int edge = 0; edge < elem->n_edges(); edge++)
    {
      // The edge's first vertex: if owned, populate coordinate vectors
      const Node & vert_node = elem->node_ref(elem->local_edge_node(edge, 0));
      const dof_id_type vert_dof = vert_node.dof_number(lagrange_sys.number(), 0, 0);

      x.set(vert_dof, vert_node(0));
      y.set(vert_dof, vert_node(1));
      z.set(vert_dof, vert_node(2));

      // The edge's second vertex: if owned, populate coordinate vectors
      const Node & wert_node = elem->node_ref(elem->local_edge_node(edge, 1));
      const dof_id_type wert_dof = wert_node.dof_number(lagrange_sys.number(), 0, 0);

      x.set(wert_dof, wert_node(0));
      y.set(wert_dof, wert_node(1));
      z.set(wert_dof, wert_node(2));

      // The edge's (middle) node: if owned, populate discrete gradient matrix
      const Node & edge_node = elem->node_ref(elem->local_edge_node(edge, 2));
      const dof_id_type edge_dof = edge_node.dof_number(nedelec_sys.number(), 0, 0);

      const Real sign = vert_node > wert_node ? 1 : -1;

      G.set(edge_dof, vert_dof,  sign);
      G.set(edge_dof, wert_dof, -sign);

    }

  // Assemble the discrete gradient matrix and the coordinate vectors
  G.close();
  x.close();
  y.close();
  z.close();

  // Compute the matrix-vector products
  G.vector_mult(Gx, x);
  G.vector_mult(Gy, y);
  G.vector_mult(Gz, z);

  // Get the preconditioner object associated with system's linear solver and hand over the matrix and vectors
  const PC & pc = dynamic_cast<PetscLinearSolver<double> &>(dynamic_cast<NewtonSolver &>(*(nedelec_sys.time_solver->diff_solver().get())).get_linear_solver()).pc();
  LibmeshPetscCall(PCHYPRESetDiscreteGradient(pc, G.mat()));
  LibmeshPetscCall(PCHYPRESetEdgeConstantVectors(pc, Gx.vec(), Gy.vec(), Gz.vec()));
}

#endif // HYPRE_AMS_H
