
```markdown
# Finite Difference Method for Heat Conduction

This project implements the finite difference method to solve the heat conduction equation in both Cartesian and cylindrical coordinates. It is designed to be run on the TCHPC Callan cluster or Chuck using MPI.

## Prerequisites

Before running the project, make sure to load the necessary modules:

```bash
module load openmpi
module load gcc
```

## Installation

To build the project, follow these steps:

<!-- 1. **Clone the repository** (if applicable):

    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ``` -->

1. **Create and navigate to the build directory**:

    ```bash
    mkdir build
    cd build
    ```

2. **Compile the project using CMake**:
    ```bash
    cmake ..
    make
    ```

## Usage

### Cartesian Coordinates

To run the simulation for Cartesian coordinates, use the following command:

```bash
mpirun -np <number-of-processes> ./main
```

### Cylindrical Coordinates

To run the simulation for cylindrical coordinates, use this command:

```bash
mpirun -np <number-of-processes> ./cylinder
```

Replace `<number-of-processes>` with the desired number of MPI processes.

## Testing

A `test.sh` script is provided to facilitate running multiple test cases and obtaining the final iteration results. To execute the script, run:

```bash
./test.sh
```

This script will handle running the simulations and saving the results for each test case.

## Visualization

Results from the simulations can be visualized using the provided MATLAB scripts:

- **`vis.m`**: For visualizing Cartesian coordinates results.
- **`cylinder.m`**: For visualizing cylindrical coordinates results.

Open MATLAB and run the appropriate script to view the results.



## Acknowledgments

- TCHPC Callan cluster and Chuck for providing computational resources.
- OpenMPI and GCC for their contributions to parallel computing and compiling.


