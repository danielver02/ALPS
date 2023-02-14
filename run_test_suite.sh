numproc=8
echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_kpar_fast_v00'
echo 'Testing: ' $input
mkdir solution
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors in generate_distribution execution of' $input
else
    echo 'ERRORS in generate_distribution execution of' $input
fi
date
cd ../
echo '-=-=-=-=-=-=-='
echo 'ALPS'
date
mpirun -np $numproc --oversubscribe ./src/ALPS $input.in 2> $input'.error' > $input'.out'
jj=$(wc -l $input'.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_kperp_v00'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors in generate_distribution execution of' $input
else
    echo 'ERRORS in generate_distribution execution of' $input
fi
date
cd ../
echo '-=-=-=-=-=-=-='
echo 'ALPS'
date
mpirun -np $numproc --oversubscribe ./src/ALPS $input.in 2> $input'.error' > $input'.out'
jj=$(wc -l $input'.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_ICW_v00'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors in generate_distribution execution of' $input
else
    echo 'ERRORS in generate_distribution execution of' $input
fi
date
cd ../
echo '-=-=-=-=-=-=-='
echo 'ALPS'
date
mpirun -np $numproc --oversubscribe ./src/ALPS $input.in 2> $input'.error' > $input'.out'
jj=$(wc -l $input'.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_map_v00'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors in generate_distribution execution of' $input
else
    echo 'ERRORS in generate_distribution execution of' $input
fi
date
cd ../
echo '-=-=-=-=-=-=-='
echo 'ALPS'
date
mpirun -np $numproc --oversubscribe ./src/ALPS $input.in 2> $input'.error' > $input'.out'
jj=$(wc -l $input'.error' | awk '{print $1}')
if [ $jj -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date
