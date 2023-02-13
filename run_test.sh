cd distribution
./generate_distribution test_kpar_v00.in
cd ../
mpirun -np 8 --oversubscribe ./src/ALPS test_kpar_fast.in
