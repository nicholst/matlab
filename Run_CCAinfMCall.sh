#!/bin/bash

Label=(
    "N00200_Ny0050_Nx0050_Nz00_Npca030_HJ0" 
    "N00200_Ny0050_Nx0050_Nz05_Npca030_HJ0" 
    "N00200_Ny0050_Nx0050_Nz05_Npca030_HJ1" 
    "N01000_Ny0200_Nx0200_Nz00_Npca050_HJ0" 
    "N01000_Ny0200_Nx0200_Nz10_Npca050_HJ0" 
    "N01000_Ny0200_Nx0200_Nz10_Npca050_HJ1" 
    "N10000_Ny2000_Nx2000_Nz00_Npca200_HJ0" 
    "N10000_Ny2000_Nx2000_Nz10_Npca200_HJ0" 
    "N10000_Ny2000_Nx2000_Nz10_Npca200_HJ1" 
)

Param=(
    "N=200;   Ny=50;   Nx=50;   Nz=0;   Npca=30;  nR=10;  nP=1000; HJ=false" 
    "N=200;   Ny=50;   Nx=50;   Nz=5;   Npca=30;  nR=10;  nP=1000; HJ=false"
    "N=200;   Ny=50;   Nx=50;   Nz=5;   Npca=30;  nR=10;  nP=1000; HJ=true"
    "N=1000;  Ny=200;  Nx=200;  Nz=0;   Npca=50;  nR=10;  nP=1000; HJ=false"
    "N=1000;  Ny=200;  Nx=200;  Nz=10;  Npca=50;  nR=10;  nP=1000; HJ=false"
    "N=1000;  Ny=200;  Nx=200;  Nz=10;  Npca=50;  nR=10;  nP=1000; HJ=true"
    "N=10000; Ny=2000; Nx=2000; Nz=0;   Npca=200; nR=1;   nP=1000; HJ=false"
    "N=10000; Ny=2000; Nx=2000; Nz=50;  Npca=200; nR=1;   nP=1000; HJ=false"
    "N=10000; Ny=2000; Nx=2000; Nz=50;  Npca=200; nR=1;   nP=1000; HJ=true"
)

BaseNm=CCAinf-

nSet=${#Param[*]}

# Changes to SGE/Slurm parameters:
#
# Small (non-Biobank) scale
#   -l h_vmem=3G      3G max mem
#   -t 1:10000        10k jobs, 10 reals. per job, total of 100k realisations 
#   nR=10             10 realisations per job
#
# Biobank scale
#   -l h_vmem=8G      8G max mem
#   -t 1:1000         1k jobs, 1 real. per job, 1k realisations... not sure how long 10k no less 100k would take
#   nR=1              1 realisations per job

for ((i=0;i<nSet;i++)) ; do
    ScriptNm=${BaseNm}${Label[i]}

    cat > $ScriptNm.sh <<EOF
#!/bin/bash
#$ -l h_rt=01:00:00
EOF

    if ((i<6)) ; then
	cat >> $ScriptNm.sh <<EOF
#$ -l h_vmem=3G
#$ -t 1:100
EOF
    else
	cat >> $ScriptNm.sh <<EOF
#$ -l h_vmem=8G
#$ -t 1:10
EOF
    fi

    cat >> $ScriptNm.sh <<EOF
#$ -cwd
. /etc/profile
module add matlab

matlab -nodisplay -nojvm -singleCompThread -r \
  "try;${Param[i]};CCAinfMC(N,Ny,Nx,Nz,Npca,nR,nP,HJ);end;quit"

EOF

    qsub -N $ScriptNm $ScriptNm.sh
done

