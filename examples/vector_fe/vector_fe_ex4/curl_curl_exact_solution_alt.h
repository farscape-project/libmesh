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

#ifndef CURL_CURL_EXACT_SOLUTION_H
#define CURL_CURL_EXACT_SOLUTION_H

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

class CurlCurlExactSolution
{
public:
  CurlCurlExactSolution() = default;
  ~CurlCurlExactSolution() = default;

  RealGradient operator()(Real x, Real y, Real z)
  {
    const Real ux = sin(k*y);
    const Real uy = sin(k*z);
    const Real uz = sin(k*x);

    return RealGradient(ux, uy, uz);
  }

  RealTensor grad(Real x, Real y, Real z)
  {
    const Real dux_dy = k*cos(k*y);
    const Real duy_dz = k*cos(k*z);
    const Real duz_dx = k*cos(k*x);

    return RealTensor(0, dux_dy, 0, 0, 0, duy_dz, duz_dx, 0, 0);
  }

  RealGradient forcing(Real x, Real y, Real z)
  {
    const Real fx = (1 + k*k)*sin(k*y);
    const Real fy = (1 + k*k)*sin(k*z);
    const Real fz = (1 + k*k)*sin(k*x);

    return RealGradient(fx, fy, fz);
  }

private:
  const Real k = .5*pi;
};

#endif // CURL_CURL_EXACT_SOLUTION_H
