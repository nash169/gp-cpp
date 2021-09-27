# Riemannian Gaussian Processes
Gaussian Processes on Manifolds

## Run
- Generate ground truth         -> `python scripts/ground_truth.py path/to/mesh` (use stl file)
- Generate samples              -> `./build/src/examples/generate_target <name>`
- FEM eigenfunction             -> `mpirun -np 1 ./build/src/examples/fem_laplace <name> <num-modes>`
- Diffusion maps eigenfunction  -> `mpirun -np 1 ./build/src/examples/diffusion_laplace <name> <num-modes>`
- Ambient solution              -> `./build/src/examples/ambient_solution <name>`
- FEM solution                  -> `./build/src/examples/fem_solution <name> <num-modes>`
- Diffusion maps solution       -> `./build/src/examples/diffusion_solution <name> <num-modes>`
- Plot eigenfunction            -> `./build/src/examples/plot_eigenfun <name> <fun-num>`
- Plot solution                 -> `./build/src/examples/plot_solution <name>`
- Plot embedding                -> `python scripts/plot_embedding.py <name> <num-modes>`

## ToDo
- Add sigma method to Gaussian Process
- Use internal Gaussian LLT decomposition of the gram matrix for GP methods
- Complete PetscVector and test it with PetscMatrix for multiple processes
- Check MPI/OPENMP how they share the cores
- MPI with shared-memory: https://stackoverflow.com/questions/64631418/passing-a-pointer-to-mpi-win-allocate-shared-through-a-wrapper
                        https://stackoverflow.com/questions/39912588/can-i-use-mpi-with-shared-memory
                        https://stackoverflow.com/questions/38592854/how-to-send-the-data-in-an-eigenmatrixxd-with-mpi