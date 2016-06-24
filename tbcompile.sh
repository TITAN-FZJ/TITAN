#!/bin/bash

# Get optional flags
while (( "$#" )); do
   case $1 in
      (uff|iff|osx|juropa|juropatest|jureca|juqueen)
        addplatform="$1"
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
        rule="$1 ${rule}"
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
      (*=*)
        addvariable="$1 ; "${addvariable}
        addvariablefilename=${addvariablefilename}"_$1"
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

cd ../../

# echo "make $rule $platform $parallel $debug $perform $filename $verbose"

# lattice="LATTICE=$(pwd | awk -F/ '{print $(NF-1)}')"
# system="SYSTEM=$(pwd | awk -F/ '{print $(NF)}')"

# Adding platform to executable filename
if [[ -z "$filename" ]] ; then
  filename=$(echo FILE=main_${addplatform}.exe)
else
  if [[ ! "$filename" =~ "$addplatform" ]] ; then
    filename="${filename:0:${#filename}-4}_"${addplatform}".exe"
  fi
fi

# Adding specific variables to mod_io.F90
if [[ ! -z "$addvariable" ]] ; then
  echo "Adding the following line to mod_io.F90:"
  echo "${addvariable}"
  nl=$'\n'
  sed -i.bak -e "/User manual additions/a\ ${nl}${addvariable}" source/mod_io.F90
  filename="${filename:0:${#filename}-4}"${addvariablefilename}".exe"
fi

# Test if "recompile"  was used before
if [[ "$rule" =~ "recompile" && ! -z `grep "recompile" last_compilation` ]] ; then
  while true; do
    read -n 1 -p "Last compilation used \"recompile\". Use again? " yn
    echo ' '
    # read -n 1 -e -p "Last compilation used \"recompile\". Use again? \n" yn
    case $yn in
      ([Yy]*)
        break
        ;;
      ([Nn]*)
        exit
        ;;
      (*)
        echo "Please answer (y)es or (n)o." ;;
    esac
  done
fi

echo "make $rule $platform $parallel $debug $perform $filename $verbose" | tee last_compilation

make $rule $platform $parallel $debug $perform $filename $verbose

# Removing added lines to mod_io.F90
if [[ ! -z "$addvariable" ]] ; then
  echo "Recovering original mod_io.F90"
  mv source/mod_io.F90.bak source/mod_io.F90
fi

# echo "error = " $?

exit