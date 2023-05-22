#!/bin/bash

echo 'Compiling RCD program...'
if [[ $(shopt nocasematch) =~ off ]] ; then
    shopt -s nocasematch
    ShoptTurnOff=true
else
    ShoptTurnOff=false
fi 

if [[ $(uname) =~ "cygwin" ]] ; then
    LAPACK=/lib/liblapack.dll.a
    BLAS=/lib/libblas.dll.a
    OUTFILE=rcd.exe
elif [[ $(uname) =~ "linux" ]] ; then
    LAPACK=-llapack
    BLAS=-lblas
    OUTFILE=rcd.out
fi


if [[ -z "$1" ]]; then
	prec=8
else
	prec=$1
fi

if [[ $3 == 'DEBUG' ]]; then
	DEBUGFLAGS='-fcheck=all -ggdb -O0 -ffpe-trap=invalid,zero,overflow'
	echo 'Will use 0th level of optimization with debugging (-fcheck=all -ggdb -O0 -ffpe-trap=invalid,zero,overflow).'
else
	DEBUGFLAGS='-O3'
	echo 'Will use 3rd level of optimization (-O3).'
fi

if [[ $2 == 'MPI' ]]; then
	FC=mpifort
	EXTRA='-lmpi'
	echo 'Will compile a parallel program. (NOT YET IMPLEMENTED)'
else
	FC=gfortran
	EXTRA=''
	echo 'Will compile a sequential program.'
fi

echo 'Do not mind the 2 (only 2!) type mismatch warnings. They will be there since I wanted to check for numerical precision errors and LAPACK does not have quadruple precision function.'

if [[ $prec == 8 ]]; then
	echo 'Will compile in double precision.'
	$FC -o $OUTFILE $DEBUGFLAGS -fallow-argument-mismatch -ffree-line-length-300 -fmax-errors=1 rcd.f95 $BLAS $LAPACK $EXTRA
elif [[ $prec == 16 ]]; then
	echo 'Will compile in quadruple precision'
	echo 'WARNING: The matrix diagonalization is still done in double precision'
	$FC -o $OUTFILE -Wall $DEBUGFLAGS -fallow-argument-mismatch -freal-8-real-16 -freal-4-real-8 -ffree-line-length-300 -fmax-errors=1 rcd.f95 $BLAS $LAPACK $EXTRA
elif [[ $prec == 4 ]]; then
	echo 'Will compile in single precision (float)'
	$FC -o $OUTFILE $DEBUGFLAGS -fallow-argument-mismatch -freal-8-real-4 -ffree-line-length-300 -fmax-errors=1 rcd.f95 $BLAS $LAPACK $EXTRA
else
	echo 'Enter 4 (single precision), 8 (double precision) or 16 (quadruple precision)'
fi


if [[ ShoptTurnOff =~ true ]] ; then
    shopt -u nocasematch
fi
