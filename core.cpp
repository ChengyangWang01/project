#ifndef CORE_CPP
#define CORE_CPP

#pragma once 
#include "heat.hpp"
#include "heat_cylinder.hpp"
#include <cstring>
#include <cmath>
#include <mpi.h>


auto inv2 = [](const double & x){ return 1.0 / (x * x); };

void 
exchange1 (
  Field_cart & field, 
  environment3d & env)
{
  int flag;
  double * sbuf, * rbuf;

  // Send / Receive ( Dim 0 )
  flag = 0;
  sbuf = &field(field.nx, 1, 1), rbuf = &field(0, 1, 1);
  MPI_Sendrecv( sbuf, 1, field.halos[0], env.neighbors_dest[0], flag, 
                rbuf, 1, field.halos[0], env.neighbors_src[0],  flag, 
                env.cart_comm, MPI_STATUS_IGNORE);

  sbuf = &field(1,1,1), rbuf = &field(field.nx+1, 1, 1);
  MPI_Sendrecv( sbuf, 1, field.halos[0], env.neighbors_src[0],  flag, 
                rbuf, 1, field.halos[0], env.neighbors_dest[0], flag, 
                env.cart_comm, MPI_STATUS_IGNORE);

  // Send / Receive ( Dim 1 )
  flag = 1;
  sbuf = &field(1, field.ny, 1), rbuf = &field(1, 0, 1);
  MPI_Sendrecv( sbuf, 1, field.halos[1], env.neighbors_dest[1], flag, 
                rbuf, 1, field.halos[1], env.neighbors_src[1],  flag, 
                env.cart_comm, MPI_STATUS_IGNORE);

  sbuf = &field(1,1,1), rbuf = &field(1, field.ny+1, 1);
  MPI_Sendrecv( sbuf, 1, field.halos[1], env.neighbors_src[1],  flag, 
                rbuf, 1, field.halos[1], env.neighbors_dest[1], flag, 
                env.cart_comm, MPI_STATUS_IGNORE);


  // // Send / Receive ( Dim 2 )
  flag = 2;
  sbuf = &field(1, 1, field.nz), rbuf = &field(1, 1, 0);
  MPI_Sendrecv( sbuf, 1, field.halos[2], env.neighbors_dest[2], flag, 
                rbuf, 1, field.halos[2], env.neighbors_src[2],  flag, 
                env.cart_comm, MPI_STATUS_IGNORE);

  sbuf = &field(1,1,1), rbuf = &field(1, 1, field.nz+1);
  MPI_Sendrecv( sbuf, 1, field.halos[2], env.neighbors_src[2],  flag, 
                rbuf, 1, field.halos[2], env.neighbors_dest[2], flag, 
                env.cart_comm, MPI_STATUS_IGNORE);

}

void exchange2(Field_cart & field, environment3d & env) {
  int flag;
  double *sbuf, *rbuf;
  MPI_Request requests[12]; 

  // Send / Receive ( Dim 0 )
  flag = 0;
  sbuf = &field(field.nx, 1, 1);
  rbuf = &field(0, 1, 1);
  MPI_Irecv(rbuf, 1, field.halos[0], env.neighbors_src[0], flag, env.cart_comm, &requests[0]);
  MPI_Isend(sbuf, 1, field.halos[0], env.neighbors_dest[0], flag, env.cart_comm, &requests[1]);

  sbuf = &field(1, 1, 1);
  rbuf = &field(field.nx + 1, 1, 1);
  MPI_Irecv(rbuf, 1, field.halos[0], env.neighbors_dest[0], flag, env.cart_comm, &requests[2]);
  MPI_Isend(sbuf, 1, field.halos[0], env.neighbors_src[0], flag, env.cart_comm, &requests[3]);

  // Send / Receive ( Dim 1 )
  flag = 1;
  sbuf = &field(1, field.ny, 1);
  rbuf = &field(1, 0, 1);
  MPI_Irecv(rbuf, 1, field.halos[1], env.neighbors_src[1], flag, env.cart_comm, &requests[4]);
  MPI_Isend(sbuf, 1, field.halos[1], env.neighbors_dest[1], flag, env.cart_comm, &requests[5]);

  sbuf = &field(1, 1, 1);
  rbuf = &field(1, field.ny + 1, 1);
  MPI_Irecv(rbuf, 1, field.halos[1], env.neighbors_dest[1], flag, env.cart_comm, &requests[6]);
  MPI_Isend(sbuf, 1, field.halos[1], env.neighbors_src[1], flag, env.cart_comm, &requests[7]);

  // Send / Receive ( Dim 2 )
  flag = 2;
  sbuf = &field(1, 1, field.nz);
  rbuf = &field(1, 1, 0);
  MPI_Irecv(rbuf, 1, field.halos[2], env.neighbors_src[2], flag, env.cart_comm, &requests[8]);
  MPI_Isend(sbuf, 1, field.halos[2], env.neighbors_dest[2], flag, env.cart_comm, &requests[9]);

  sbuf = &field(1, 1, 1);
  rbuf = &field(1, 1, field.nz + 1);
  MPI_Irecv(rbuf, 1, field.halos[2], env.neighbors_dest[2], flag, env.cart_comm, &requests[10]);
  MPI_Isend(sbuf, 1, field.halos[2], env.neighbors_src[2], flag, env.cart_comm, &requests[11]);

  MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);
}



double
update(
        Field_cart & curr, 
  const Field_cart & prev,
  const double & a, 
  const double & dt
)
{
  double diff {0.0};

  auto inv_dx2 {inv2(prev.dx)}, inv_dy2 {inv2(prev.dy)}, inv_dz2 {inv2(prev.dz)};

  for (std::size_t i = 1; i < curr.nx+1; ++i)
  {
    for (std::size_t j = 1; j < curr.ny+1; ++j)
    {
      if ( i + prev.sx - prev.nx_full / 2 == 0 && j + prev.sy - prev.ny_full / 2 == 0)
      continue;
      
      for (std::size_t k = 1; k < curr.nz+1; ++k)
      {
    curr(i, j, k) = prev(i, j, k) + a * dt * 
    (
( prev(i + 1, j, k) - 2.0 * prev(i, j, k) + prev(i - 1, j, k) ) * inv_dx2 +
( prev(i, j + 1, k) - 2.0 * prev(i, j, k) + prev(i, j - 1, k) ) * inv_dy2 + 
( prev(i, j, k + 1) - 2.0 * prev(i, j, k) + prev(i, j, k - 1) ) * inv_dz2
    );
    diff += std::pow(curr(i,j,k) - prev(i,j,k), 2);
      }
    }
  }

  return diff;
}

double update(
  Field_cylinder & curr,
  Field_cylinder & prev,
  const double & a,
  const double & dt
)
{
  double diff {0.0};

  auto inv_drho2 {inv2(prev.drho)}, inv_dtheta2 {inv2(prev.dtheta)}, inv_dz2 {inv2(prev.dz)};
  auto inv_drho  {1.0 / prev.drho};

  for (std::size_t i = 1; i < curr.nrho+1; ++i)
  {
    double rho { prev.radius * (i + prev.srho - 1) / (prev.nrho_full-2)};
    for (std::size_t j = 1; j < curr.ntheta+1; ++j)
    {
      for (std::size_t k = 1; k < curr.nz+1; ++k)
      {
        curr(i,j,k) = prev(i,j,k) + a * dt * 
        (
  ( prev(i + 1, j, k) - 2.0 * prev(i, j, k) + prev(i - 1, j, k) ) * inv_drho2 
+ ( prev(i + 1, j, k)                       - prev(i - 1, j, k) ) / ( 2 * rho * prev.drho)
+ ( prev(i, j + 1, k) - 2.0 * prev(i, j, k) + prev(i, j - 1, k) ) * inv_dtheta2 * inv2(rho) * inv2(rho)
+ ( prev(i, j, k + 1) - 2.0 * prev(i, j, k) + prev(i, j, k - 1) ) * inv_dz2   
        );
diff += std::pow(curr(i,j,k) - prev(i,j,k), 2);
      }
    }
  }
  return diff;
}




double update2(
    Field_cart & curr, 
    Field_cart & prev,
    const double & a, 
    const double & dt)
{
    double diff {0.0};
    auto inv_dx2 {inv2(prev.dx)}, inv_dy2 {inv2(prev.dy)}, inv_dz2 {inv2(prev.dz)};                     
                              
    const double omega{1.2};// Set the relaxation factor omega, usually choose 1 < omega < 2.

    for (std::size_t i = 1; i < curr.nx + 1; ++i) {
        for (std::size_t j = 1; j < curr.ny + 1; ++j) {
            for (std::size_t k = 1; k < curr.nz + 1; ++k) {
                double x_minus, x_plus, y_minus, y_plus, z_minus, z_plus;

                // X-dimension
                if (i > 1) {
                    x_minus = curr(i - 1, j, k);
                } else {
                    x_minus = prev(i - 1, j, k);
                }
                x_plus = prev(i + 1, j, k);

                // Y-dimension
                if (j > 1) {
                    y_minus = curr(i, j - 1, k);
                } else {
                    y_minus = prev(i, j - 1, k);
                }
                y_plus = prev(i, j + 1, k);

                // Z-dimension
                if (k > 1) {
                    z_minus = curr(i, j, k - 1);
                } else {
                    z_minus = prev(i, j, k - 1);
                }
                z_plus = prev(i, j, k + 1);

                double gauss_seidel_value = prev(i, j, k) + a * dt * (
                    (x_plus - 2.0 * prev(i, j, k) + x_minus) * inv_dx2 +
                    (y_plus - 2.0 * prev(i, j, k) + y_minus) * inv_dy2 +
                    (z_plus - 2.0 * prev(i, j, k) + z_minus) * inv_dz2);

                // sor
                double old_value = curr(i, j, k);
                curr(i, j, k) = (1 - omega) * old_value + omega * gauss_seidel_value;
                diff += std::pow(curr(i, j, k) - old_value, 2);
            }
        }
    }

    return diff;
}




void 
gather(
  Field_cart       & subField,
  Matrix3D<double> & globField,
  environment3d    & env
)
{

  if (env.size == 1) { // num process = 1; just do a  swap;
    globField.swap(subField.temperature);
  } else {

  MPI_Datatype sub_type;

  int pid, i, j, k;
  int xidx {1}, yidx {1}, zidx {1};

  int sxs[env.size], sys[env.size], szs[env.size];
  int nxs[env.size], nys[env.size], nzs[env.size];

  int sxyz_cpy[3] = {subField.sx, subField.sy, subField.sz};
  int nx {subField.nx}, ny {subField.ny}, nz {subField.nz};

  if (subField.sx == 1) { -- sxyz_cpy[0]; -- xidx; ++ nx; }
  if (subField.sy == 1) { -- sxyz_cpy[1]; -- yidx; ++ ny; }
  if (subField.sz == 1) { -- sxyz_cpy[2]; -- zidx; ++ nz; }

  if (subField.ex == globField.nx-2) ++ nx;
  if (subField.ey == globField.ny-2) ++ ny;
  if (subField.ez == globField.nz-2) ++ nz;

  MPI_Gather(&sxyz_cpy[0], 1, MPI_INT, sxs, 1, MPI_INT, 0, env.cart_comm);
  MPI_Gather(&sxyz_cpy[1], 1, MPI_INT, sys, 1, MPI_INT, 0, env.cart_comm);
  MPI_Gather(&sxyz_cpy[2], 1, MPI_INT, szs, 1, MPI_INT, 0, env.cart_comm);

  MPI_Gather(&nx, 1, MPI_INT, nxs, 1, MPI_INT, 0, env.cart_comm);
  MPI_Gather(&ny, 1, MPI_INT, nys, 1, MPI_INT, 0, env.cart_comm);
  MPI_Gather(&nz, 1, MPI_INT, nzs, 1, MPI_INT, 0, env.cart_comm);

  ///
  if (env.rank != 0) {
    int array_size[]    = {subField.nx+2, subField.ny+2, subField.nz+2};
    int array_subsize[] = {nx, ny, nz}; int array_start[] = {0, 0, 0};
    
    MPI_Type_create_subarray( 3,  array_size, 
                                  array_subsize, 
                                  array_start,
                                  MPI_ORDER_C, MPI_DOUBLE, &sub_type);  

    MPI_Type_commit(&sub_type);
    MPI_Send(&subField(xidx, yidx, zidx), 1, sub_type, 0, env.rank, env.cart_comm);
    MPI_Type_free(&sub_type);
  }

  if (env.rank == 0) {
    for (pid = 0; pid < env.size; ++pid) {
      if (pid == 0) {
        for (i = sxyz_cpy[0]; i <= subField.ex+1; ++i) 
        {
          for (j = sxyz_cpy[1]; j <= subField.ey; ++j) 
          {
  memcpy( &globField(i,j,sxyz_cpy[2]),
          &subField(i,j,sxyz_cpy[2]), 
          nzs[pid]*sizeof(double));
          }
        }
      }

      if (pid != 0) {
        int array_subsize[] = {nxs[pid], nys[pid], nzs[pid]}, array_start[] = {0, 0, 0};
        int array_size[] = {globField.nx, globField.ny, globField.nz};

        MPI_Type_create_subarray(3, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_DOUBLE, &sub_type);

        MPI_Type_commit(&sub_type);
        MPI_Recv(&globField(sxs[pid], sys[pid], szs[pid]), 1, sub_type, pid, pid, env.cart_comm, MPI_STATUS_IGNORE);
        MPI_Type_free(&sub_type);
      }
    }
  }


  }


}

#endif // end define CORE_CPP