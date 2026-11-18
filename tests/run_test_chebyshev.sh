nnodes=10
input='test_chebyshev'
cd ..
mkdir -p solution
cd interpolation
./interpolation test_interp.in
cd ../distribution
./generate_distribution $input'_dist.in'
cd ../
mpirun -np $nnodes --oversubscribe ./src/ALPS 'tests/'$input'.in'

input='test_chebyshev'
cd ..
mkdir -p solution
cd interpolation
./interpolation test_interp.in
cd ../distribution
./generate_distribution $input'_dist.in'
cd ../
mpirun -np $nnodes --oversubscribe ./src/ALPS 'tests/'$input'.in'
