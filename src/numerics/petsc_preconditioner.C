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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

// Local Includes
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"

namespace libMesh
{

template <typename T>
PetscPreconditioner<T>::PetscPreconditioner (const libMesh::Parallel::Communicator & comm_in) :
  Preconditioner<T>(comm_in)
{}



template <typename T>
void PetscPreconditioner<T>::apply(const NumericVector<T> & x, NumericVector<T> & y)
{
  PetscVector<T> & x_pvec = cast_ref<PetscVector<T> &>(const_cast<NumericVector<T> &>(x));
  PetscVector<T> & y_pvec = cast_ref<PetscVector<T> &>(const_cast<NumericVector<T> &>(y));

  Vec x_vec = x_pvec.vec();
  Vec y_vec = y_pvec.vec();

  PetscErrorCode ierr = PCApply(_pc, x_vec, y_vec);
  LIBMESH_CHKERR(ierr);
}




template <typename T>
void PetscPreconditioner<T>::init ()
{
  libmesh_error_msg_if(!this->_matrix, "ERROR: No matrix set for PetscPreconditioner, but init() called");

  // Clear the preconditioner in case it has been created in the past
  if (!this->_is_initialized)
    {
      // Should probably use PCReset(), but it's not working at the moment so we'll destroy instead
      if (_pc)
        _pc.destroy();

      PetscErrorCode ierr = PCCreate(this->comm().get(), _pc.get());
      LIBMESH_CHKERR(ierr);

      auto pmatrix = cast_ptr<PetscMatrixBase<T> *>(this->_matrix);
      _mat = pmatrix->mat();
    }

  PetscErrorCode ierr = PCSetOperators(_pc, _mat, _mat);
  LIBMESH_CHKERR(ierr);

  // Set the PCType.  Note: this used to be done *before* the call to
  // PCSetOperators(), and only when !_is_initialized, but
  // 1.) Some preconditioners (those employing sub-preconditioners,
  // for example) have to call PCSetUp(), and can only do this after
  // the operators have been set.
  // 2.) It should be safe to call set_petsc_preconditioner_type()
  // multiple times.
  set_petsc_preconditioner_type(this->_preconditioner_type, *_pc);

  this->_is_initialized = true;
}



template <typename T>
void PetscPreconditioner<T>::clear()
{
  // Calls custom deleter
  _pc.destroy();
}



template <typename T>
PC PetscPreconditioner<T>::pc()
{
  return _pc;
}



template <typename T>
void PetscPreconditioner<T>::set_petsc_preconditioner_type (const PreconditionerType & preconditioner_type, PC & pc)
{
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // get the communicator from the PETSc object
  Parallel::communicator comm;
  ierr = PetscObjectGetComm((PetscObject)pc, & comm);
  if (ierr != LIBMESH_PETSC_SUCCESS)
    libmesh_error_msg("Error retrieving communicator");
  Parallel::Communicator communicator(comm);

  switch (preconditioner_type)
    {
    case IDENTITY_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCNONE));
      CHKERRABORT(comm,ierr);
      break;

    case CHOLESKY_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCCHOLESKY));
      CHKERRABORT(comm,ierr);
      break;

    case ICC_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCICC));
      CHKERRABORT(comm,ierr);
      break;

    case ILU_PRECOND:
      {
        // In serial, just set the ILU preconditioner type
        if (communicator.size())
          {
            ierr = PCSetType (pc, const_cast<KSPType>(PCILU));
            CHKERRABORT(comm,ierr);
          }
        else
          {
            // But PETSc has no truly parallel ILU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign ILU sub-preconditioners.
            ierr = PCSetType (pc, const_cast<KSPType>(PCBJACOBI));
            CHKERRABORT(comm,ierr);

            // Set ILU as the sub preconditioner type
            set_petsc_subpreconditioner_type(PCILU, pc);
          }
        break;
      }

    case LU_PRECOND:
      {
        // In serial, just set the LU preconditioner type
        if (communicator.size())
          {
            ierr = PCSetType (pc, const_cast<KSPType>(PCLU));
            CHKERRABORT(comm,ierr);
          }
        else
          {
            // But PETSc has no truly parallel LU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign LU sub-preconditioners.
            ierr = PCSetType (pc, const_cast<KSPType>(PCBJACOBI));
            CHKERRABORT(comm,ierr);

            // Set ILU as the sub preconditioner type
            set_petsc_subpreconditioner_type(PCLU, pc);
          }
        break;
      }

    case ASM_PRECOND:
      {
        // In parallel, I think ASM uses ILU by default as the sub-preconditioner...
        // I tried setting a different sub-preconditioner here, but apparently the matrix
        // is not in the correct state (at this point) to call PCSetUp().
        ierr = PCSetType (pc, const_cast<KSPType>(PCASM));
        CHKERRABORT(comm,ierr);
        break;
      }

    case JACOBI_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCJACOBI));
      CHKERRABORT(comm,ierr);
      break;

    case BLOCK_JACOBI_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCBJACOBI));
      CHKERRABORT(comm,ierr);
      break;

    case SOR_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCSOR));
      CHKERRABORT(comm,ierr);
      break;

    case EISENSTAT_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCEISENSTAT));
      CHKERRABORT(comm,ierr);
      break;

    case AMG_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCHYPRE));
      CHKERRABORT(comm,ierr);
      break;

    case SVD_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCSVD));
      CHKERRABORT(comm,ierr);
      break;

    case USER_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCMAT));
      CHKERRABORT(comm,ierr);
      break;

    case SHELL_PRECOND:
      ierr = PCSetType (pc, const_cast<KSPType>(PCSHELL));
      CHKERRABORT(comm,ierr);
      break;

    default:
      libMesh::err << "ERROR:  Unsupported PETSC Preconditioner: "
                   << preconditioner_type       << std::endl
                   << "Continuing with PETSC defaults" << std::endl;
    }

  // Set additional options if we are doing AMG and
  // HYPRE is available
#ifdef LIBMESH_HAVE_PETSC_HYPRE
  if (preconditioner_type == AMG_PRECOND)
    {
      ierr = PCHYPRESetType(pc, "boomeramg");
      CHKERRABORT(comm,ierr);
    }
#endif

  // Let the commandline override stuff
  ierr = PCSetFromOptions(pc);
  CHKERRABORT(comm,ierr);
}


template <typename T>
void PetscPreconditioner<T>::set_petsc_subpreconditioner_type(const PCType type, PC & pc)
{
  // For catching PETSc error return codes
  PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;

  // get the communicator from the PETSc object
  Parallel::communicator comm;
  ierr = PetscObjectGetComm((PetscObject)pc, & comm);
  if (ierr != LIBMESH_PETSC_SUCCESS)
    libmesh_error_msg("Error retrieving communicator");
  Parallel::Communicator communicator(comm);

  // All docs say must call KSPSetUp or PCSetUp before calling PCBJacobiGetSubKSP.
  // You must call PCSetUp after the preconditioner operators have been set, otherwise you get the:
  //
  // "Object is in wrong state!"
  // "Matrix must be set first."
  //
  // error messages...
  ierr = PCSetUp(pc);
  CHKERRABORT(comm,ierr);

  // To store array of local KSP contexts on this processor
  KSP * subksps;

  // the number of blocks on this processor
  PetscInt n_local;

  // The global number of the first block on this processor.
  // This is not used, so we just pass null instead.
  // int first_local;

  // Fill array of local KSP contexts
  ierr = PCBJacobiGetSubKSP(pc, &n_local, LIBMESH_PETSC_NULLPTR,
                            &subksps);
  CHKERRABORT(comm,ierr);

  // Loop over sub-ksp objects, set ILU preconditioner
  for (PetscInt i=0; i<n_local; ++i)
    {
      // Get pointer to sub KSP object's PC
      PC subpc;
      ierr = KSPGetPC(subksps[i], &subpc);
      CHKERRABORT(comm,ierr);

      // Set requested type on the sub PC
      ierr = PCSetType(subpc, type);
      CHKERRABORT(comm,ierr);
    }
}


#ifdef LIBMESH_HAVE_PETSC_HYPRE
template <typename T>
void PetscPreconditioner<T>::set_petsc_hypre_aux_data(PC & pc, System * sys)
{
  // Get the hypre preconditioner we are using
  const char * hypre_type = nullptr;
  LibmeshPetscCall(PCSetFromOptions(pc));
  LibmeshPetscCall(PCHYPREGetType(pc, &hypre_type));

  // If not ams or not a 1st order Nédélec system, we quit as we do not support anything else at the moment
  if (!hypre_type || std::string(hypre_type) != "ams" || !sys || sys->n_vars() != 1 || sys->variable(0).type() != FEType(1, NEDELEC_ONE))
    return;

  // The mesh we build the preconditioner for
  const MeshBase & mesh = sys->get_equation_systems().get_mesh();

  // Dummy Lagrange system defined over the same mesh so we can enumerate the vertices
  System & lagrange_sys = sys->get_equation_systems().add_system<System>("dummy");
  lagrange_sys.add_variable("dummy");
  lagrange_sys.reinit_mesh();

  // Global (i.e. total) and local (i.e. to this processor) number of edges and vertices
  const dof_id_type n_glb_edges = sys->n_dofs();
  const dof_id_type n_loc_edges = sys->n_local_dofs();
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

      if (vert_node.processor_id() == global_processor_id())
      {
        x.set(vert_dof, vert_node(0));
        y.set(vert_dof, vert_node(1));
        z.set(vert_dof, vert_node(2));
      }

      // The edge's second vertex: if owned, populate coordinate vectors
      const Node & wert_node = elem->node_ref(elem->local_edge_node(edge, 1));
      const dof_id_type wert_dof = wert_node.dof_number(lagrange_sys.number(), 0, 0);

      if (wert_node.processor_id() == global_processor_id())
      {
        x.set(wert_dof, wert_node(0));
        y.set(wert_dof, wert_node(1));
        z.set(wert_dof, wert_node(2));
      }

      // The edge's (middle) node: if owned, populate discrete gradient matrix
      const Node & edge_node = elem->node_ref(elem->local_edge_node(edge, 2));
      const dof_id_type edge_dof = edge_node.dof_number(sys->number(), 0, 0);

      if (edge_node.processor_id() == global_processor_id())
      {
        const Real sign = vert_node > wert_node ? 1 : -1;

        G.set(edge_dof, vert_dof,  sign);
        G.set(edge_dof, wert_dof, -sign);
      }
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
  LibmeshPetscCall(PCHYPRESetDiscreteGradient(pc, G.mat()));
  LibmeshPetscCall(PCHYPRESetEdgeConstantVectors(pc, Gx.vec(), Gy.vec(), Gz.vec()));
}
#endif


//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscPreconditioner<Number>;

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
