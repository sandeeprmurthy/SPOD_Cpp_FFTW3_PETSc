module restore
module load phdf5/1.14.3 
module load eigen/3.4.0
module load fftw3/3.3.10

export PETSC_DIR=/home1/03119/srmurth2/lib/petsc-3.21.2
export PETSC_ARCH=arch-linux-c-debug
export SLEPC_DIR=/home1/03119/srmurth2/lib/slepc-3.21.1

export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/apps/intel24/impi21/phdf5/1.14.3/x86_64/lib:/opt/apps/intel24/impi21/fftw3/3.3.10/lib:$LD_LIBRARY_PATH
export CPATH=/opt/apps/intel24/impi21/phdf5/1.14.3/x86_64/include:/opt/apps/intel24/impi21/fftw3/3.3.10/include:$CPATH

rm ./spod_program
make
./spod_program
#ibrun ./spod_program

#rm ./read_data
#make
#./read_data
