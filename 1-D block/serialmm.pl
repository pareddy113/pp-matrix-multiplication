#!usr/bin/perl
system('module load tools/mpich2-1.5-gcc');
system('rm serialmm.out serialmm.err');
system('mpicc -o serialmm serialmm.c');
system('qsub serialmm.pbs');

