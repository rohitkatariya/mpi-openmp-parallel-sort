#!/bin/sh
### Set the job name (for your reference)
#PBS -N rohitcol730
### Set the project name, your department code by default
#PBS -P cse
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M csz208844@iitd.ac.in
####
###PBS -l select=16:ncpus=3:mpiprocs=1
#PBS -l select=3:ncpus=4:mpiprocs=1:ompthreads=2
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=00:10:00
#PBS -P col730.csz208844  
#PBS -o jobOutput.txt
#PBS -l place=scatter
###PBS -l software=OpenMPI
# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired. 
#echo "==============================="
echo $PBS_JOBID
echo "PBS_NTASKS:".$PBS_NTASKS
#cat $PBS_NODEFILE
#echo "==============================="
echo "omp num threads echo ".$OMP_NUM_THREADS
cd $PBS_O_WORKDIR
#job 
#time -p ls 
module load compiler/gcc/6.5/openmpi/4.0.2
#mpirun -n 4 a.out 
#time -p mpirun -np $PBS_NTASKS a.out /home/cse/phd/csz208844/file_transfer/chhavi/out_file> $PBS_JOBID
#export OMP_NUM_THREADS=4
time -p mpirun a.out input_dir/inputfile100000> $PBS_JOBID

#NOTE
# The job line is an example : users need to change it to suit their applications
# The PBS select statement picks n nodes each having m free processors
# OpenMPI needs more options such as $PBS_NODEFILE
