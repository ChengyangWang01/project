#ifndef HEAT_HPP
#define HEAT_HPP

#include "base.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>


/// @brief Basic MPI cartesian topology environment structure
struct environment3d
{
  int size, rank;
  const int dim {3};
  int dims[3] = {0,0,0}, periods[3], coords[3] = {0,0,0};
  int neighbors_src[3], neighbors_dest[3]; 
  MPI_Comm cart_comm;

  environment3d(int periods[3]) : periods {periods[3]}
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Dims_create(size, dim, dims);
    MPI_Cart_create(MPI_COMM_WORLD, dim, dims, periods, 0, &cart_comm);

    for (int i = 0; i < dim; ++i)
      MPI_Cart_shift(cart_comm, i, 1, &neighbors_src[i], &neighbors_dest[i]);
    
    MPI_Cart_coords(cart_comm, rank, dim, coords);

    // ====================== A few information for DEBUGGING ====================== //
    std::cout << "Constructing 3D mpi cartesian environment \n";
    std::cout << " Rank : " 
    << std::fixed << std::setw(3) << rank << "/" 
    << std::fixed << std::setw(3) << size;

    std::cout << " has Neighbors: \n";
    for (int i = 0; i < dim; ++i)
      std::cout << "\t\tDirection " << i 
      << " - [" 
      << std::fixed << std::setw(3) << neighbors_src[i] 
      << ", " 
      << std::fixed << std::setw(3) << neighbors_dest[i] << "]." << std::endl;
  }

}; // end of mpi environment3d

struct Field_cart
{
  typedef int                   size_type;
  typedef const int             const_size_type;
  typedef double                value_type;
  typedef Matrix3D<double>      array_type;

  array_type    temperature;

  size_type     nx, ny, nz;
  size_type     sx, sy, sz, ex, ey, ez;
  size_type     nx_full, ny_full, nz_full;

  value_type    dx {0.01}, dy {0.01}, dz {0.01};

  value_type    radius { 1.0 / 3.0 };

  MPI_Datatype  halos[3];

  Field_cart(const_size_type, const_size_type, const_size_type, environment3d &);

  void inits(environment3d &);

  double& operator()(const_size_type, const_size_type, const_size_type);
  const double& operator()(const_size_type, const_size_type, const_size_type) const;


  ~Field_cart() { std::cout << "Calling Destructor of Cartesian Field." << std::endl; }
}; // end of struct Field Cartesian



inline 
Field_cart::Field_cart(
  const_size_type nx_in, 
  const_size_type ny_in, 
  const_size_type nz_in, 
  environment3d & env)
{
  nx_full = nx_in; ny_full = ny_in;  nz_full = nz_in; 

  auto decomp = [](const int & n, const int & size, const int & rank, 
    size_type &n_loc, size_type &start, size_type &end)
  {
    n_loc = n / (( size_type ) size) ;
    start = n_loc * rank + 1;
    size_type deficit { n % (( size_type ) size) } ;

    start = start + ( ( rank < deficit ) ? rank : deficit );
    end = start + n_loc - 1;
    if ( end > n || rank == size-1 ) end = n;
    if ( rank < deficit ) ++n_loc ;
    return 0;
  };

  decomp(nx_full-2, env.dims[0], env.coords[0], nx, sx, ex);
  decomp(ny_full-2, env.dims[1], env.coords[1], ny, sy, ey);
  decomp(nz_full-2, env.dims[2], env.coords[2], nz, sz, ez);

  temperature = Matrix3D<double>(nx+2, ny+2, nz+2);
  
  // ================================= Create Halos ================================= //
  int array_sizes[]     = {nx+2, ny+2, nz+2}, 
      array_starts[]    = {0,0,0}, 
      array_subsizes[]  = {nx, ny, nz};

  MPI_Type_vector(ny, nz, nz+2, MPI_DOUBLE, &halos[0]);                // Dim 0 

  array_subsizes[1] = 1;
  MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts,
                            MPI_ORDER_C, MPI_DOUBLE, &halos[1]);           // Dim 1

  array_subsizes[1] = ny; 
  array_subsizes[2] = 1; 
  MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts,
                            MPI_ORDER_C, MPI_DOUBLE, &halos[2]);           // Dim 2

  for (auto & halo : halos) { MPI_Type_commit(&halo); }



  // ====================== A few information for DEBUGGING ====================== //

  // std::cout << "Rank : " << env.rank << "/" << env.size << "\n"
  //           << "Constructing subarray with shape : "
  //           << "[" << nx <<", " << ny << ", " << nz << "]"
  //           << " of full shape "
  //           << "[" << nx_full <<", " << ny_full << ", " << nz_full << "]" << "\n"
  //           << "\t[" << sx <<", " << sy << ", " << sz << "]" << " - "
  //           << "[" << ex <<", " << ey << ", " << ez << "]"
  //           << std::endl;

}

inline
double& 
Field_cart::operator()(const_size_type i, const_size_type j, const_size_type k)
{ return temperature(i,j,k); }


inline
const double& 
Field_cart::operator()(const_size_type i, const_size_type j, const_size_type k)
const { return temperature(i,j,k); }



inline 
void 
Field_cart::inits(environment3d & env)
{
  auto r { nx_full * radius };
  
  // Fill inner parts
  for (size_type i = 0; i < nx+2; ++i)
  {
    for (size_type j = 0; j < ny+2; ++j)
    {
      auto dx { i + sx - nx_full / 2 };
      auto dy { j + sy - ny_full / 2 };
      for (size_type k = 0; k < nz+2; ++k)
      {
        if (dx * dx + dy * dy < r * r) 
        {
          temperature(i,j,k) = 0.0;
        } else if (dx * dx + dy * dy == 0)
        {
          temperature(i,j,k) = 0.0;
        } else {
          temperature(i,j,k) = 65.0;
        } 
      }
    }
  }
}

#endif // end define HEAT_HPP