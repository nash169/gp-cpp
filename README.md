# Riemannian Gaussian Processes
Gaussian Processes on Manifolds

## Run
- Generate ground truth         -> `python scripts/ground_truth path/to/mesh` (use stl file)
- Generate samples              -> `./build/srx/examples/generate_target` (set mesh name within the exe)
- FEM eigenfunction             -> `mpirun -np 1 ./build/srx/examples/fem_laplace` (set mesh name within the exe)
- Diffusion maps eigenfunction  -> `mpirun -np 1 ./build/srx/examples/diffusion_laplace` (set mesh name within the exe)
- Ambient solution              -> `./build/srx/examples/ambient_solution` (set mesh name within the exe)
- FEM solution                  -> `./build/srx/examples/fem_solution` (set mesh name within the exe)
- Diffusion maps solution       -> `./build/srx/examples/diffusion_solution` (set mesh name within the exe)