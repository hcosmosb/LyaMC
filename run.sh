OMP_NUM_THREADS=2
export OMP_NUM_THREADS
gfortran -fwrapv -o main.exe rng.f90 main.f90 initial_setup.f90 run.f90 LyART.f90 CalcT_sub.f90 Hyojeong.f90 exit_condition.f90 PO_plane.f90 -fopenmp && ./main.exe 1
