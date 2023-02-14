mkdir solution
cd interpolation
./interpolation test_interp.in
cd ../distribution
./generate_distribution test_kpar_fast_v00_dist.in
cd ../
mpirun -np 8 --oversubscribe ./src/ALPS test_kpar_fast_v00.in
