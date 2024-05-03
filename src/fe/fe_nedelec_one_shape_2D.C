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


// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// An excellent discussion of Nedelec shape functions is given in
// https://www.dealii.org/reports/nedelec/nedelec.pdf
template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  const short sign = elem->point(i) > elem->point((i+1) % elem->n_vertices()) ? 1 : -1;

  const Real xi  = p(0);
  const Real eta = p(1);

  switch (total_order)
    {
      // linear Nedelec (first kind) shape functions
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              libmesh_assert_less (i, 4);

              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              switch(i)
                {
                case 0:
                  return sign * RealGradient( -0.25*(1.0-eta), 0.0 );
                case 1:
                  return sign * RealGradient( 0.0, -0.25*(1.0+xi) );
                case 2:
                  return sign * RealGradient( 0.25*(1.0+eta), 0.0 );
                case 3:
                  return sign * RealGradient( 0.0, -0.25*(xi-1.0) );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }

              return RealGradient();
            }

          case TRI6:
          case TRI7:
            {
              libmesh_assert_less (i, 3);

              switch(i)
                {
                case 0:
                  return sign * RealGradient( -1.0+eta, -xi );
                case 1:
                  return sign * RealGradient( eta, -xi );
                case 2:
                  return sign * RealGradient( eta, -xi+1.0 );

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }


    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD9:
            {
              libmesh_assert_less (i, 12);

              const Real xi  = p(0);
              const Real eta = p(1);

              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              Real x = 0.5 * (xi + 1.0);
              Real y = 0.5 * (eta + 1.0);

              switch(i)
                {
                case 0:
                  {
                    if (elem->point(0) > elem->point(1))
                        return RealGradient(-9.0*x*y*y+12*x*y-3*x+6*y*y-8*y+2, 0.0);
                    else
                        return RealGradient( 9.0*x*y*y-12*x*y+3*x-6*y*y+8*y-2, 0.0);
                  }
                case 1:
                  {
                    if (elem->point(0) > elem->point(1))
                      return RealGradient( 9.0*x*y*y-12*x*y+3*x-3*y*y+4*y-1, 0.0);
                    else
                      return RealGradient(-9.0*x*y*y+12*x*y-3*x+3*y*y-4*y+1, 0.0);
                  }

                case 2:
                  {
                    if (elem->point(1) > elem->point(2))
                        return RealGradient(0.0,  x*(-9*x*y+6*x+6*y-4));
                    else
                        return RealGradient(0.0,  x*( 9*x*y-6*x-6*y+4));
                  }
                case 3:
                  {
                    if (elem->point(1) > elem->point(2))
                      return RealGradient(0.0, x*( 9*x*y-3*x-6*y+2));
                    else
                      return RealGradient(0.0, x*(-9*x*y+3*x+6*y-2));
                  }

                case 4:
                  {
                    if (elem->point(2) < elem->point(3))
                      return RealGradient(y*(-9.0*x*y+6*x+6*y-4), 0.0);
                    else
                      return RealGradient(y*( 9.0*x*y-6*x-6*y+4), 0.0);
                  }
                case 5:
                  {
                    if (elem->point(2) < elem->point(3))
                      return RealGradient(y*( 9.0*x*y-6*x-3*y+2), 0.0);
                    else
                      return RealGradient(y*(-9.0*x*y+6*x+3*y-2), 0.0);
                  }

                case 6:
                  {
                    if (elem->point(3) < elem->point(0))
                      return RealGradient(0.0, -9*x*x*y+6*x*x+12*x*y-8*x-3*y+2);
                    else
                      return RealGradient(0.0,  9*x*x*y-6*x*x-12*x*y+8*x+3*y-2);
                  }

                case 7:
                  {
                    if (elem->point(3) < elem->point(0))
                      return RealGradient(0.0,  9*x*x*y-3*x*x-12*x*y+4*x+3*y-1);
                    else
                      return RealGradient(0.0, -9*x*x*y+3*x*x+12*x*y-4*x-3*y+1);
                  }

                case 8:
                    return RealGradient(0.0,  3*x*(3*x*y-2*x-3*y+2));
                
                case 9:
                    return RealGradient(3*y*(-3*x*y+3*x+2*y-2), 0.0);
                
                case 10:
                    return RealGradient(3*y*(3*x*y-3*x-y+1), 0.0);
                
                case 11:
                    return RealGradient(0.0,  3*x*(-3*x*y+x+3*y-1));

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }

              return RealGradient();
            }

          case TRI6:
          case TRI7:
            {
              const Real xi  = p(0);
              const Real eta = p(1);

              libmesh_assert_less (i, 3);

              switch(i)
                {
                case 0:
                  {
                    if (elem->point(0) > elem->point(1))
                      return RealGradient( -1.0+eta, -xi);
                    else
                      return RealGradient( 1.0-eta, xi);
                  }
                case 1:
                  {
                    if (elem->point(1) > elem->point(2))
                      return RealGradient( eta, -xi);
                    else
                      return RealGradient( -eta, xi);
                  }

                case 2:
                  {
                    if (elem->point(2) > elem->point(0))
                      return RealGradient( eta, -xi+1.0);
                    else
                      return RealGradient( -eta, xi-1.0);
                  }

                default:
                  libmesh_error_msg("Invalid i = " << i);
                }
            }
            
          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }

      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(elem, order, i, p, add_p_level);
  libmesh_not_implemented();
#endif
}




template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const FEType fet,
                                      const Elem * elem,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape(elem, fet.order, i, p, add_p_level);
}

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const Elem * elem,
                                            const Order order,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);
  libmesh_assert_less (j, 2);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());
  const short sign = elem->point(i) > elem->point((i+1) % elem->n_vertices()) ? 1 : -1;

  switch (total_order)
    {
      // linear Nedelec (first kind) shape function first derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              libmesh_assert_less (i, 4);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                      case 2:
                        return RealGradient();
                      case 1:
                      case 3:
                        return sign * RealGradient( 0.0, -0.25 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 1:
                      case 3:
                        return RealGradient();
                      case 0:
                      case 2:
                        return sign * RealGradient( 0.25 );

                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }

          case TRI6:
          case TRI7:
            {
              libmesh_assert_less (i, 3);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return sign * RealGradient( 0.0, -1.0 );

                  // d()/deta
                case 1:
                  return sign * RealGradient( 1.0 );

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }
    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD9:
            {
              libmesh_assert_less (i, 12);
              
              const Real xi  = p(0);
              const Real eta = p(1);

              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              Real x = 0.5 * (xi + 1.0);
              Real y = 0.5 * (eta + 1.0);

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                       case 0:
                         {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.5*(-9*y*y+12*y-3), 0.0);
                            else
                              return RealGradient(0.5*( 9*y*y-12*y+3), 0.0);
                          }
                      case 1:
                        {
                          if (elem->point(0) > elem->point(1))
                            return RealGradient(0.5*( 9.0*y*y-12*y+3), 0.0);
                          else
                            return RealGradient(0.5*(-9.0*y*y+12*y-3), 0.0);
                        }
                      case 2:
                       {
                          if (elem->point(1) > elem->point(2))
                            return RealGradient(0.0, 0.5*(-18*x*y+12*x+6*y-4));
                          else
                            return RealGradient(0.0, 0.5*( 18*x*y-12*x-6*y+4));
                        }
                      case 3:
                        {
                          if (elem->point(1) > elem->point(2))
                            return RealGradient(0.0, 0.5*( 18*x*y-6*x-6*y+2));
                          else
                            return RealGradient(0.0, 0.5*(-18*x*y+6*x+6*y-2));
                        }
                    case 4:
                        {
                          if (elem->point(2) < elem->point(3))
                            return RealGradient(0.5*(-9.0*y*y+6*y), 0.0 );
                          else
                            return RealGradient(0.5*( 9.0*y*y-6*y), 0.0 );
                        }
                    case 5:
                        {
                          if (elem->point(2) < elem->point(3))
                            return RealGradient(0.5*( 9.0*y*y-6*y), 0.0);
                          else
                            return RealGradient(0.5*(-9.0*y*y+6*y), 0.0 );
                        }
                    case 6:
                       {
                        if (elem->point(3) < elem->point(0))
                          return RealGradient(0.0, 0.5*(-18*x*y+12*x+12*y-8));
                        else
                          return RealGradient(0.0, 0.5*( 18*x*y-12*x-12*y+8));
                       }
                    case 7:
                       {
                        if (elem->point(3) < elem->point(0))
                          return RealGradient(0.0,  0.5*( 18*x*y-6*x-12*y+4));
                        else
                          return RealGradient(0.0,  0.5*(-18*x*y+6*x+12*y-4));
                       }
                    case 8:
                      return RealGradient(0.0,  1.5*(6*x*y-4*x-3*y+2));

                    case 9:
                      return RealGradient(1.5*y*(-3*y+3), 0.0);
    
                    case 10:
                      return RealGradient(1.5*y*(3*y-3), 0.0);
                
                    case 11:
                      return RealGradient(0.0,  1.5*(-6*x*y+2*x+3*y-1));
                      
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                       case 0:
                         {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.5*(-18.0*x*y+12*x+12*y-8), 0.0);
                            else
                              return RealGradient(0.5*( 18.0*x*y-12*x-12*y+8), 0.0);
                          }
                        case 1:
                          {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.5*( 18.0*x*y-12*x-6*y+4), 0.0);
                            else
                              return RealGradient(0.5*(-18.0*x*y+12*x+6*y-4), 0.0);
                          }
                        case 2:
                          {
                            if (elem->point(1) > elem->point(2))
                              return RealGradient(0.0, 0.5*x*(-9*x+6));
                            else
                              return RealGradient(0.0, 0.5*x*( 9*x-6));
                          }
                        case 3:
                          {
                            if (elem->point(1) > elem->point(2))
                              return RealGradient(0.0, 0.5*x*( 9*x-6));
                            else
                              return RealGradient(0.0, 0.5*x*(-9*x+6));
                          }
                        case 4:
                          {
                            if (elem->point(2) < elem->point(3))
                              return RealGradient(0.5*(-18.0*x*y+6*x+12*y-4), 0.0 );
                            else
                              return RealGradient(0.5*( 18.0*x*y-6*x-12*y+4), 0.0 );
                          }
                        case 5:
                          {
                            if (elem->point(2) < elem->point(3))
                              return RealGradient(0.5*( 18.0*x*y-6*x-6*y+2), 0.0);
                            else
                              return RealGradient(0.5*(-18.0*x*y+6*x+6*y-2), 0.0);
                          }
                        case 6:
                          {
                           if (elem->point(3) < elem->point(0))
                            return RealGradient(0.0, 0.5*(-9*x*x+12*x-3));
                          else
                            return RealGradient(0.0, 0.5*( 9*x*x-12*x+3));
                         }
                        case 7:
                          {
                            if (elem->point(3) < elem->point(0))
                              return RealGradient(0.0,  0.5*( 9*x*x-12*x+3));
                            else
                              return RealGradient(0.0,  0.5*(-9*x*x+12*x-3));
                          }
                        case 8:
                          return RealGradient(0.0,  1.5*x*(3*x-3));

                        case 9:
                          return RealGradient(1.5*(-6*x*y+3*x+4*y-2), 0.0);
  
                       case 10:
                          return RealGradient(1.5*(6*x*y-3*x-2*y+1), 0.0);
                
                        case 11:
                          return RealGradient(0.0, 1.5*x*(-3*x+3));
                     
                        default:
                          libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }

              return RealGradient();
            }

          case TRI6:
          case TRI7:
            {
              libmesh_assert_less (i, 3);

              // Account for edge flipping
              Real f = 1.0;

              switch(i)
                {
                case 0:
                  {
                    if (elem->point(0) > elem->point(1))
                      f = -1.0;
                    break;
                  }
                case 1:
                  {
                    if (elem->point(1) > elem->point(2))
                      f = -1.0;
                    break;
                  }
                case 2:
                  {
                    if (elem->point(2) > elem->point(0))
                      f = -1.0;
                    break;
                  }
                default:
                  libmesh_error_msg("Invalid i = " << i);
                }

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    return RealGradient( 0.0, f*1.0);
                  }
                  // d()/deta
                case 1:
                  {
                    return RealGradient( f*(-1.0) );
                  }
                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
            }


          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));
          }
      }
      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);
    }
#else // LIBMESH_DIM > 1
  libmesh_ignore(elem, order, i, j, add_p_level);
  libmesh_not_implemented();
#endif
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const ElemType,
                                            const Order,
                                            const unsigned int,
                                            const unsigned int,
                                            const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const FEType fet,
                                            const Elem * elem,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}





#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const Elem * elem,
                                                   const Order order,
                                                   const unsigned int libmesh_dbg_var(i),
                                                   const unsigned int libmesh_dbg_var(j),
                                                   const Point & p,
                                                   const bool add_p_level)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  libmesh_assert_less (j, 3);

  const Order total_order = static_cast<Order>(order + add_p_level * elem->p_level());

  switch (total_order)
    {
      // linear Nedelec (first kind) shape function second derivatives
    case FIRST:
      {
        switch (elem->type())
          {
          case QUAD8:
          case QUAD9:
            {
              libmesh_assert_less (i, 4);
              // All second derivatives for linear quads are zero.
              return RealGradient();
            }

          case TRI6:
          case TRI7:
            {
              libmesh_assert_less (i, 3);
              // All second derivatives for linear triangles are zero.
              return RealGradient();
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));

          } // end switch (type)
      } // end case FIRST


    case SECOND:
      {
        switch (elem->type())
          {
          case QUAD9:
            {
              libmesh_assert_less (i, 12);
              
              const Real xi  = p(0);
              const Real eta = p(1);

              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
              libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
              libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

              Real x = 0.5 * (xi + 1.0);
              Real y = 0.5 * (eta + 1.0);

              switch (j)
                {
                  //d^2 () / dxi^2
                case 0:
                  {
                    switch(i)
                      {
                       case 0:
                         {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.0, 0.0);
                            else
                              return RealGradient(0.0, 0.0);
                          }
                      case 1:
                        {
                          if (elem->point(0) > elem->point(1))
                            return RealGradient(0.0, 0.0);
                          else
                            return RealGradient(0.0, 0.0);
                        }
                      case 2:
                        {
                          if (elem->point(1) > elem->point(2))
                            return RealGradient(0.0, 0.25*(-18*y+12));
                          else
                            return RealGradient(0.0, 0.25*( 18*y-12));
                        }
                      case 3:
                        {
                          if (elem->point(1) > elem->point(2))
                            return RealGradient(0.0, 0.25*( 18*y-6));
                          else
                            return RealGradient(0.0, 0.25*(-18*y+6));
                        }
                    case 4:
                        {
                          if (elem->point(2) < elem->point(3))
                            return RealGradient(0.0, 0.0 );
                          else
                            return RealGradient(0.0, 0.0 );
                        }
                    case 5:
                        {
                          if (elem->point(2) < elem->point(3))
                            return RealGradient(0.0, 0.0);
                          else
                            return RealGradient(0.0, 0.0 );
                        }
                    case 6:
                       {
                        if (elem->point(3) < elem->point(0))
                          return RealGradient(0.0, 0.25*(-18*y+12));
                        else
                          return RealGradient(0.0, 0.25*( 18*y-12));
                       }
                    case 7:
                       {
                        if (elem->point(3) < elem->point(0))
                          return RealGradient(0.0,  0.25*( 18*y-6));
                        else
                          return RealGradient(0.0,  0.25*(-18*y+6));
                       }
                    case 8:
                      return RealGradient(0.0,  0.75*(6*y-4));

                    case 9:
                      return RealGradient(0.0, 0.0);
    
                    case 10:
                      return RealGradient(0.0, 0.0);
                
                    case 11:
                      return RealGradient(0.0,  0.75*(-6*y+2));
                      
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 0

                  // ^2()/dxi deta
                case 1:
                  {
                    switch(i)
                      {
                       case 0:
                         {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.25*(-18*y+12), 0.0);
                            else
                              return RealGradient(0.25*( 18*y-12), 0.0);
                          }
                      case 1:
                        {
                          if (elem->point(0) > elem->point(1))
                            return RealGradient(0.25*( 18.0*y-12), 0.0);
                          else
                            return RealGradient(0.25*(-18.0*y+12), 0.0);
                        }
                      case 2:
                       {
                          if (elem->point(1) > elem->point(2))
                            return RealGradient(0.0, 0.25*(-18*x+6));
                          else
                            return RealGradient(0.0, 0.25*( 18*x-6));
                        }
                      case 3:
                        {
                          if (elem->point(1) > elem->point(2))
                            return RealGradient(0.0, 0.25*( 18*x-6));
                          else
                            return RealGradient(0.0, 0.25*(-18*x+6));
                        }
                    case 4:
                        {
                          if (elem->point(2) < elem->point(3))
                            return RealGradient(0.25*(-18.0*y+6), 0.0);
                          else
                            return RealGradient(0.25*( 18.0*y-6), 0.0);
                        }
                    case 5:
                        {
                          if (elem->point(2) < elem->point(3))
                            return RealGradient(0.25*( 18.0*y-6), 0.0);
                          else
                            return RealGradient(0.25*(-18.0*y+6), 0.0);
                        }
                    case 6:
                       {
                        if (elem->point(3) < elem->point(0))
                          return RealGradient(0.0, 0.25*(-18*x+12));
                        else
                          return RealGradient(0.0, 0.25*( 18*x-12));
                       }
                    case 7:
                       {
                        if (elem->point(3) < elem->point(0))
                          return RealGradient(0.0,  0.25*( 18*x-12));
                        else
                          return RealGradient(0.0,  0.25*(-18*x+12));
                       }
                    case 8:
                      return RealGradient(0.0,  0.75*(6*x-3));

                    case 9:
                      return RealGradient(0.75*(-6*y), 0.0);
    
                    case 10:
                      return RealGradient(0.75*(6*y), 0.0);
                
                    case 11:
                      return RealGradient(0.0, 0.75*(-6*x+3));
                      
                      default:
                        libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 1

                    //d^2 () / deta^2
                case 2:
                  {
                    switch(i)
                      {
                       case 0:
                         {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.25*(-18.0*x+12), 0.0);
                            else
                              return RealGradient(0.25*( 18.0*x-12), 0.0);
                          }
                        case 1:
                          {
                            if (elem->point(0) > elem->point(1))
                              return RealGradient(0.25*( 18.0*x-6), 0.0);
                            else
                              return RealGradient(0.25*(-18.0*x+6), 0.0);
                          }
                        case 2:
                          {
                            if (elem->point(1) > elem->point(2))
                              return RealGradient(0.0,  0.0);
                            else
                              return RealGradient(0.0,  0.0);
                          }
                        case 3:
                          {
                            if (elem->point(1) > elem->point(2))
                              return RealGradient(0.0, 0.0);
                            else
                              return RealGradient(0.0, 0.0);
                          }
                        case 4:
                          {
                            if (elem->point(2) < elem->point(3))
                              return RealGradient(0.25*(-18.0*x+12), 0.0 );
                            else
                              return RealGradient(0.25*( 18.0*x-12), 0.0 );
                          }
                        case 5:
                          {
                            if (elem->point(2) < elem->point(3))
                              return RealGradient(0.25*( 18.0*x-6), 0.0);
                            else
                              return RealGradient(0.25*(-18.0*x+6), 0.0);
                          }
                        case 6:
                          {
                           if (elem->point(3) < elem->point(0))
                            return RealGradient(0.0, 0.0);
                          else
                            return RealGradient(0.0, 0.0);
                         }
                        case 7:
                          {
                            if (elem->point(3) < elem->point(0))
                              return RealGradient(0.0, 0.0);
                            else
                              return RealGradient(0.0, 0.0);
                          }
                        case 8:
                          return RealGradient(0.0, 0.0);

                        case 9:
                          return RealGradient(0.75*(-6*x+4), 0.0);
  
                       case 10:
                          return RealGradient(0.75*(6*x-2), 0.0);
                
                        case 11:
                          return RealGradient(0.0, 0.0);
                     
                        default:
                          libmesh_error_msg("Invalid i = " << i);
                      }
                  } // j = 2

                default:
                  libmesh_error_msg("Invalid j = " << j);
                }
              return RealGradient();
            }

          case TRI6:
          case TRI7:
            {
              libmesh_assert_less (i, 3);
              // All second derivatives for linear triangles are zero.
              return RealGradient();
            }

          default:
            libmesh_error_msg("ERROR: Unsupported 2D element type!: " << Utility::enum_to_string(elem->type()));

          } // end switch (type)
      } // end case FIRST




      // unsupported order
    default:
      libmesh_error_msg("ERROR: Unsupported 2D FE order!: " << total_order);

    } // end switch (order)

#else // LIBMESH_DIM > 1
  libmesh_assert(true || i || j);
  libmesh_ignore(elem, order, add_p_level);
  libmesh_not_implemented();
#endif
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const ElemType,
                                                   const Order,
                                                   const unsigned int,
                                                   const unsigned int,
                                                   const Point &)
{
  libmesh_error_msg("Nedelec elements require the element type \nbecause edge orientation is needed.");
  return RealGradient();
}


template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const FEType fet,
                                                   const Elem * elem,
                                                   const unsigned int i,
                                                   const unsigned int j,
                                                   const Point & p,
                                                   const bool add_p_level)
{
  return FE<2,NEDELEC_ONE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}



#endif

} // namespace libMesh
