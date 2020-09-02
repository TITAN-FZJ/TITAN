#This script should be executed with two arguments
#such arguments should be paths to the files we want to compare
#for example
#$ ./compare_output_files.sh folder1/output/selfconsistency folder2/output/selfconsistency

# for the moment it does not distinguish between 0.0 and -0.0


file1="$1"
file2="$2"

printf "\nI am comparing the files $1 and $2\n"

###############################################
printf "Fermi energy\n"

dummy="Self-consistent ground state"

diff <(grep "$dummy" -A 1 $file1) <(grep "$dummy" -A 1 $file2) > /dev/null

if [ $? -eq 0 ]; then
    printf ">>>The Fermi energies are identical\n"
else
    diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2)
fi
###############################################
printf "\nCharge density\n"

dummy="Charge density"

diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2) > /dev/null

if [ $? -eq 0 ]; then
    printf ">>>The charge densities are identical\n"
else
    diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2)
fi
###############################################
printf "\nGap parameters\n"

dummy="Averages of the norms per orbital type"

diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2) > /dev/null

if [ $? -eq 0 ]; then
    printf ">>>The gap parameters are identical\n"
else
    printf ">>>These results are different, it can be problem of -0.0 != 0.0\n"
    #diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2)
fi
###############################################
printf "\nMagnetization components\n"

dummy="Magnetization components"

diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2) > /dev/null

if [ $? -eq 0 ]; then
    printf ">>>The magnetization components are identical\n"
else
    printf ">>>These results are different, it can be problem of -0.0 != 0.0\n"
    #diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2)
fi
###############################################
printf "\nOrbital components\n"

dummy="Orbital components in global frame"

diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2) > /dev/null

if [ $? -eq 0 ]; then
    printf ">>>The orbital components are identical\n"
else
    printf ">>>These results are different, it can be problem of -0.0 != 0.0\n"
    #diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2)
fi
###############################################
printf "\nTotal\n"

dummy="Total"

diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2) > /dev/null

if [ $? -eq 0 ]; then
    printf ">>>The total results are identical\n"
else
    diff <(grep "$dummy" -A 11 $file1) <(grep "$dummy" -A 11 $file2)
fi
###############################################
printf "\n"
dummy="Time after self-consistency"
grep "$dummy" $file1
grep "$dummy" $file2
printf "\n"
dummy="Finished on"
grep "$dummy" $file1
grep "$dummy" $file2
###############################################
