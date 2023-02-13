cd distribution
./generate_distribution.e test_kpar_v00.in
cd ../
mpirun -np 8 ./ALPS.e test_kpar_fast.in
