#!/bin/bash

# Get optional flags
while (( "$#" )); do
   case $1 in
      (uff|iff|osx|juropa|juropatest|jureca|juqueen)
        platform="PLATFORM=$1"
        shift
        ;;
      (omp)
        parallel="PARALLEL=$1"
        shift
        ;;
      (debug)
        debug="DEBUG=$1"
        shift
        ;;
      (scalasca)
        perform="PERFORM=$1"
        shift
        ;;
      (cleandep|cleandebug|cleanobj|cleanmod|cleanexe|clean|cleanall|recompile)
        rule="$1"
        shift
        ;;
      (*.exe)
        if [[ -z "${1:0:${#1}-4}" ]] ; then
          echo "Empty filename!"
          exit 1
        fi
        filename=$(echo FILE=$1)
        shift
        ;;
      (verbose)
        verbose="--debug=$1"
        shift
        ;;
      (*)
        echo "Illegal option: $1"
        exit 1
   esac
done

# echo "make $platform $parallel $debug $filename"

# ln ../../makefile ./build 2>/dev/null
# if [ -z $rule ] ; then
#   ln ../../source/*.F90 ../source/*.F90 ./source/*.F90 ../../f90_mod_deps.py ./build 2>/dev/null
# fi
# cd ./build

lattice="LATTICE=$(pwd | awk -F/ '{print $(NF-1)}')"
system="SYSTEM=$(pwd | awk -F/ '{print $(NF)}')"

cd ../../
echo "make $lattice $system $rule $platform $parallel $debug $perform $filename $verbose"
make $lattice $system $rule $platform $parallel $debug $perform $filename $verbose
# cd -