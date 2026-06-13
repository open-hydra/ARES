#!/bin/bash -
#===============================================================================
#
#          FILE: MOSE.sh
#
#         USAGE: run "./MOSE.sh [options]" in the current shell
#
#   DESCRIPTION: A script to compile and run MOSE
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     compile accordingly with the presets file"
  echo " solve                 -->     run MOSE"
  echo " kill                  -->     kill the process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> OpenMP threads"
  echo " -m | --mpi <n>        -->     launch solver with <n> MPI processes"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=./../../
MASTER=$MASTERDIR/bin/ARES
LOCAL=./bin/ARES

# Default Options
BG=0
NTHREADS=1
NPROCS=1

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case $1 in

    -b | --background)
        #echo " -> Background"
        BG=1
        shift
        ;;
    -p | --parallel)
        #echo " -> OpenMP($2)"
        NTHREADS=$2
        shift 2
        ;;
    -m | --mpi)
        #echo " -> MPI($2)"
        NPROCS=$2
        shift 2
        ;;
    -h | --help)
        print_usage
        shift
        ;;
    * )
      break
      ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi
   
if [[ $1 == solve ]]; then
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS

  # Binding: ogni rank MPI occupa NTHREADS core contigui
  export OMP_PROC_BIND=close
  export OMP_PLACES=cores

  # Check executable
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi

  # Costruisci il comando base
  if [[ $NPROCS -gt 1 ]]; then
    CMD="mpirun -np $NPROCS --map-by socket:PE=$NTHREADS $LOCAL"
  else
    CMD="$LOCAL"
  fi

  # Run the solver
  if [[ $BG == 0 ]]; then
    $CMD
  else
    $CMD 2>errors_file >logfile &
    echo $! > .ID
  fi
fi

if [[ $1 == kill ]]; then
read PID < .ID && kill $PID
fi
