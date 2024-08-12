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
    const Real ux = sin(kx*x + ky*y + kz*z);
    const Real uy = sin(kx*x + ky*y + kz*z);
    const Real uz = sin(kx*x + ky*y + kz*z);

    return RealGradient(ux, uy, uz);
  }

  RealTensor grad(Real x, Real y, Real z)
  {
    const Real du_dx = kx*cos(kx*x + ky*y + kz*z);
    const Real du_dy = ky*cos(kx*x + ky*y + kz*z);
    const Real du_dz = kz*cos(kx*x + ky*y + kz*z);

    return RealTensor(du_dx, du_dy, du_dz, du_dx, du_dy, du_dz, du_dx, du_dy, du_dz);
  }

  RealGradient forcing(Real x, Real y, Real z)
  {
    const Real fx = (1 + ky*ky + kz*kz - kx*(ky + kz))*sin(kx*x + ky*y + kz*z);
    const Real fy = (1 + kx*kx + kz*kz - ky*(kx + kz))*sin(kx*x + ky*y + kz*z);
    const Real fz = (1 + kx*kx + ky*ky - kz*(kx + ky))*sin(kx*x + ky*y + kz*z);

    return RealGradient(fx, fy, fz);
  }

private:
  const Real kx = .3*pi, ky = .4*pi, kz = .5*pi;
};

#endif // CURL_CURL_EXACT_SOLUTION_H
