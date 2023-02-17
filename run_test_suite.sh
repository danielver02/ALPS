echo '-=-=-=-=-=-=-='
echo 'Interpolate Test'
echo '-=-=-=-=-=-=-='
date
cd interpolation
./interpolation test_interp.in 2> 'test_interp.error' > 'test_interp.out'
ll0=$(wc -l 'test_interp.error' | awk '{print $1}')
if [ $ll0 -eq 0 ]
then
    echo 'No Errors in Interpolation Test'
else
    echo 'ERRORS in Interpolation Test'
fi
date
cd ../

numproc=12
echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_kpar_fast'
echo 'Testing: ' $input
mkdir -p solution
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj0=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj0 -eq 0 ]
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
kk0=$(wc -l $input'.error' | awk '{print $1}')
if [ $kk0 -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date


echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_kperp'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj1=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj1 -eq 0 ]
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
kk1=$(wc -l $input'.error' | awk '{print $1}')
if [ $kk1 -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_ICW'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj2=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj2 -eq 0 ]
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
kk2=$(wc -l $input'.error' | awk '{print $1}')
if [ $kk2 -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_map'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj3=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj3 -eq 0 ]
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
kk3=$(wc -l $input'.error' | awk '{print $1}')
if [ $kk3 -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_relativistic'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj4=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj4 -eq 0 ]
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
kk4=$(wc -l $input'.error' | awk '{print $1}')
if [ $kk4 -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

echo '-=-=-=-=-=-=-='
echo ''
echo '-=-=-=-=-=-=-='
input='test_bimax'
echo 'Testing: ' $input
echo '-=-=-=-=-=-=-='
echo 'Generate Distribution'
date
cd distribution
./generate_distribution $input'_dist.in' 2> $input'_dist.error' > $input'_dist.out'
jj5=$(wc -l $input'_dist.error' | awk '{print $1}')
if [ $jj5 -eq 0 ]
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
kk5=$(wc -l $input'.error' | awk '{print $1}')
if [ $kk5 -eq 0 ]
then
    echo 'No Errors for ALPS Execution of ' $input
else
    echo "ERRORS in ALPS Execution of" $input
fi
date

if [ $jj0 -eq 0 ] && [ $kk0 -eq 0 ] && [ $jj1 -eq 0 ] && [ $kk1 -eq 0 ] && [ $jj2 -eq 0 ] && [ $kk2 -eq 0 ] && [ $jj3 -eq 0 ] && [ $kk3 -eq 0 ] && [ $jj4 -eq 0 ] && [ $kk4 -eq 0 ] && [ $jj5 -eq 0 ] && [ $kk5 -eq 0 ]
then
    echo "No errors in ALPS Execution of Test Suite. Congratulations."
else
    echo "ERRORS in ALPS Execution of Test Suite. Investigate..."
fi
