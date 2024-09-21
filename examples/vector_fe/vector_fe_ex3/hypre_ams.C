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

#include "hypre_ams.h"

HypreAMS::HypreAMS(FEMSystem & sys)
  : nedelec_sys(sys),
    lagrange_sys(sys.get_equation_systems().add_system<System>("dummy")),
    mesh(sys.get_mesh()),
    pc(dynamic_cast<PetscLinearSolver<double> &>(dynamic_cast<NewtonSolver &>(*(sys.time_solver->diff_solver().get())).get_linear_solver()).pc())
{
  Size();
  Create();
  Build();
  Mult();
  Set();
};

void HypreAMS::Size()
{
  lagrange_sys.add_variable("dummy");
  lagrange_sys.get_equation_systems().init();

  n_glb_edges = nedelec_sys.n_dofs();
  n_loc_edges = nedelec_sys.n_local_dofs();
  n_glb_vertices = lagrange_sys.n_dofs();
  n_loc_vertices = lagrange_sys.n_local_dofs();
}

void HypreAMS::Create()
{
  ierr = MatCreateAIJ(mesh.comm().get(), n_loc_edges, n_loc_vertices, n_glb_edges, n_glb_vertices, 2, NULL, 2, NULL, &G);

  ierr = VecCreateMPI(mesh.comm().get(), n_loc_vertices, n_glb_vertices, &x);
  ierr = VecCreateMPI(mesh.comm().get(), n_loc_vertices, n_glb_vertices, &y);
  ierr = VecCreateMPI(mesh.comm().get(), n_loc_vertices, n_glb_vertices, &z);

  ierr = VecCreateMPI(mesh.comm().get(), n_loc_edges, n_glb_edges, &Gx);
  ierr = VecCreateMPI(mesh.comm().get(), n_loc_edges, n_glb_edges, &Gy);
  ierr = VecCreateMPI(mesh.comm().get(), n_loc_edges, n_glb_edges, &Gz);
}

void HypreAMS::Build()
{
  for (const auto & elem : mesh.active_local_element_ptr_range())
    for(unsigned int edge = 0; edge < elem->n_edges(); edge++)
    {
      const Node & vert_node = elem->node_ref(elem->local_edge_node(edge, 0));
      const dof_id_type vert_dof = vert_node.dof_number(lagrange_sys.number(), 0, 0);

      if (vert_node.processor_id() == global_processor_id())
      {
        ierr = VecSetValue(x, PetscInt(vert_dof), PetscScalar(vert_node(0)), INSERT_VALUES);
        ierr = VecSetValue(y, PetscInt(vert_dof), PetscScalar(vert_node(1)), INSERT_VALUES);
        ierr = VecSetValue(z, PetscInt(vert_dof), PetscScalar(vert_node(2)), INSERT_VALUES);
      }

      const Node & wert_node = elem->node_ref(elem->local_edge_node(edge, 1));
      const dof_id_type wert_dof = wert_node.dof_number(lagrange_sys.number(), 0, 0);

      if (wert_node.processor_id() == global_processor_id())
      {
        ierr = VecSetValue(x, PetscInt(wert_dof), PetscScalar(wert_node(0)), INSERT_VALUES);
        ierr = VecSetValue(y, PetscInt(wert_dof), PetscScalar(wert_node(1)), INSERT_VALUES);
        ierr = VecSetValue(z, PetscInt(wert_dof), PetscScalar(wert_node(2)), INSERT_VALUES);
      }

      const Node & edge_node = elem->node_ref(elem->local_edge_node(edge, 2));
      const dof_id_type edge_dof = edge_node.dof_number(nedelec_sys.number(), 0, 0);

      if (edge_node.processor_id() == global_processor_id())
      {
        const double sign = vert_node > wert_node ? 1 : -1;
        const PetscScalar data[] = {sign, -sign};

        const PetscInt row_index[] = {PetscInt(edge_dof)};
        const PetscInt col_index[] = {PetscInt(vert_dof), PetscInt(wert_dof)};

        ierr = MatSetValues(G, 1, row_index, 2, col_index, data, INSERT_VALUES);
      }
    }

  ierr = VecAssemblyBegin(x);
  ierr = VecAssemblyBegin(y);
  ierr = VecAssemblyBegin(z);
  ierr = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);

  ierr = VecAssemblyEnd(x);
  ierr = VecAssemblyEnd(y);
  ierr = VecAssemblyEnd(z);
  ierr = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
}

void HypreAMS::Mult()
{
  ierr = MatMult(G, x, Gx);
  ierr = MatMult(G, y, Gy);
  ierr = MatMult(G, z, Gz);
}

void HypreAMS::Set()
{
  ierr = PCHYPRESetDiscreteGradient(pc, G);
  ierr = PCHYPRESetEdgeConstantVectors(pc, Gx, Gy, Gz);
}
