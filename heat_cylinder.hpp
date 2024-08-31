#ifndef HEAT_CYLINDER_HPP
#define HEAT_CYLINDER_HPP

#pragma once
#include "base.hpp"
#include "heat.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>

struct Field_cylinder : Field_cart
{
  Field_cylinder() = default;
  Field_cylinder(const_size_type, const_size_type, const_size_type, environment3d &);

  void inits(environment3d &);

  ~Field_cylinder();


  size_type & nrho      {this->nx}      , & ntheta      {this->ny};
  size_type & nrho_full {this->nx_full} , & ntheta_full {this->ny_full};
  size_type & srho      {this->sx}      , & stheta      {this->sy};

  value_type & drho     {this->dx}      , & dtheta      {this->dy};
};

inline Field_cylinder::Field_cylinder(const_size_type nrho_in, const_size_type ntheta_in, const_size_type nz_in, environment3d & env)
: Field_cart(nrho_in, ntheta_in, nz_in, env) 
{
  drho    = this->radius / (nrho_full - 2.0);
  dtheta  = 2.0 * std::numbers::pi / (ntheta_full - 2.0);
  dz      = 1.0 / (nz_full - 2.0);
}

inline Field_cylinder::~Field_cylinder()
{
  std::cout << "Calling Destructor of Cylinder Field. (Based on Cartesian)." << std::endl;
}


inline void Field_cylinder::inits(environment3d & env)
{
  for (size_type i = 0; i < this->nx+2; ++i)
  {
    for (size_type j = 0; j < this->ny+2; ++j)
    {
      for (size_type k = 0; k < this->nz+2; ++k)
      {
        if ( i + this->sx == this->nx_full )
        {
          temperature(i,j,k) = 65.0;
        } else {
          temperature(i,j,k) = 0.0;
        }
      }
    }
  }
}


#endif // end define HEAT_CYLINDER_HPP