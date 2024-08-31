

#include "heat_cylinder.hpp"
#include "io.hpp"
#include "core.cpp"
#include <iostream>


int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);
  int periods[3] {0, 1, 0};
  auto env {environment3d(periods)};

  constexpr int nmax {20};
  constexpr std::size_t nr {50}, ntheta {65}, nz {80};
  constexpr std::size_t nr_full {nr+2}, ntheta_full {ntheta+2}, nz_full {nz+2};
  double time {0.0};
  int i {0};

  Matrix3D<double> mat (nr_full, ntheta_full, nz_full);

  Field_cylinder curr (nr_full, ntheta_full, nz_full, env), prev (nr_full, ntheta_full, nz_full, env);
  
  curr.inits(env);
  prev.inits(env);

  // Diffusion constant 
  double ldiff {0.0}, gdiff {0.0};
  constexpr double a = 1, tol {1E-1};

  // Largest stable time step 
  auto inv_drho2 {inv2(curr.drho)}, inv_dtheta2 {inv2(curr.dtheta)}, inv_dz2 {inv2(curr.dz)};
  double dt =  1.0 / (
    inv_drho2 + inv_drho2 * inv_dtheta2 + inv_dz2
  ) / a * 0.5;

  auto t0 { MPI_Wtime() };
  for (i = 1; i < 10000000; ++i)
  {
    ldiff = update(curr, prev, a, dt);
    exchange1(curr, env);
 

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
    mat.write_to_file("cylinder.txt");
  }



  MPI_Finalize();
  return 0;
}