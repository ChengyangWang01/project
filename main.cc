#include "heat.hpp"
#include "heat_cylinder.hpp"
#include "io.hpp"
#include "core.cpp"
#include <iostream>

int main ( int argc, char ** argv )
{
  MPI_Init(&argc, &argv);

  constexpr int nmax {1'00000000};
  constexpr std::size_t nx {301}, ny {90}, nz {80}; 
  constexpr std::size_t nx_full {nx+2}, ny_full {ny+2}, nz_full {nz+2};
  double time {0.0};
  int i {0};

  Matrix3D<double> mat (nx_full, ny_full, nz_full);
  // std::cout << Mat << std::endl;

  int periods[3] {0, 0, 0};
  auto env {environment3d(periods)};

  Field_cart curr(nx_full, ny_full, nz_full, env), prev(nx_full, ny_full, nz_full, env) ;

  curr.inits(env);
  prev.inits(env);

  // Diffusion constant 
  double ldiff {0.0}, gdiff {0.0};
  constexpr double a = 1, tol {1E-1};

  // Largest stable time step 
  auto inv_dx2 {inv2(curr.dx)}, inv_dy2 {inv2(curr.dy)}, inv_dz2 {inv2(curr.dz)};
  auto dt = (1 / (inv_dx2 + inv_dy2 + inv_dz2)) / a * 0.5;


  auto t0 { MPI_Wtime() };
  for (i = 1; i <= nmax; ++i)
  {
    ldiff = update(curr, prev, a, dt);
    exchange2(curr, env);

    MPI_Allreduce(&ldiff, &gdiff, 1, MPI_DOUBLE, MPI_SUM, env.cart_comm);
    if (env.rank == 0 && i % 100 == 0) {
      std::cout << i << " " << gdiff << std::endl;
    }

    if (gdiff <= tol) { break; }
    curr.temperature.swap(prev.temperature);
  }
  auto t1 { MPI_Wtime() };


  double loc_time = t1 - t0;
  MPI_Reduce(&loc_time, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  gather(curr, mat, env);

  if (env.rank == 0) 
  {
    std::cout << "Converge : " << i
              << "Time: "      << time 
              << std::endl; 
    mat.write_to_file("cart.txt");
  }
  MPI_Finalize();
  return 0;
}