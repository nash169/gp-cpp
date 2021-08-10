# Riemannian Gaussian Processes
Gaussian Processes on Manifolds

## Run
- Generate ground truth         -> `python scripts/ground_truth.py path/to/mesh` (use stl file)
- Generate samples              -> `./build/src/examples/generate_target` (set mesh name within the exe)
- FEM eigenfunction             -> `mpirun -np 1 ./build/src/examples/fem_laplace` (set mesh name within the exe)
- Diffusion maps eigenfunction  -> `mpirun -np 1 ./build/src/examples/diffusion_laplace` (set mesh name within the exe)
- Ambient solution              -> `./build/src/examples/ambient_solution` (set mesh name within the exe)
- FEM solution                  -> `./build/src/examples/fem_solution` (set mesh name within the exe)
- Diffusion maps solution       -> `./build/src/examples/diffusion_solution` (set mesh name within the exe)
- Plot eigenfunction            -> `./build/src/examples/plot_eigenfun` (set mesh name within the exe)
- Plot solution                 -> `./build/src/examples/plot_solution` (set mesh name within the exe)
- Plot embedding                -> `python scripts/plot_embedding.py` (set mesh name within the exe)

## ToDo
- Add sigma method to Gaussian Process
- Use internal Gaussian LLT decomposition of the gram matrix for GP methods
- Complete PetscVector and test it with PetscMatrix for multiple processes
- Check MPI/OPENMP how they share the cores
- MPI with shared-memory: https://stackoverflow.com/questions/64631418/passing-a-pointer-to-mpi-win-allocate-shared-through-a-wrapper
                        https://stackoverflow.com/questions/39912588/can-i-use-mpi-with-shared-memory
                        https://stackoverflow.com/questions/38592854/how-to-send-the-data-in-an-eigenmatrixxd-with-mpi