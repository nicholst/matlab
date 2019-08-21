#!/bin/bash

Label=(
    "N00200_Nx0050_Ny0050_Npca030_Nz00_HJ0" 
    "N00200_Nx0050_Ny0050_Npca030_Nz05_HJ0" 
    "N00200_Nx0050_Ny0050_Npca030_Nz05_HJ1" 
    "N01000_Nx0200_Ny0200_Npca050_Nz00_HJ0" 
    "N01000_Nx0200_Ny0200_Npca050_Nz10_HJ0" 
    "N01000_Nx0200_Ny0200_Npca050_Nz10_HJ1" 
    "N10000_Nx2000_Ny2000_Npca200_Nz00_HJ0" 
    "N10000_Nx2000_Ny2000_Npca200_Nz10_HJ0" 
    "N10000_Nx2000_Ny2000_Npca200_Nz10_HJ1" 
)

Param=(
    "N=200;   Nx=50;   Ny=50;   Npca=30;  Nz=0;   nR=10;  nP=1000; HJ=false" 
    "N=200;   Nx=50;   Ny=50;   Npca=30;  Nz=5;   nR=10;  nP=1000; HJ=false"
    "N=200;   Nx=50;   Ny=50;   Npca=30;  Nz=5;   nR=10;  nP=1000; HJ=true"
    "N=1000;  Nx=200;  Ny=200;  Npca=50;  Nz=0;   nR=10;  nP=1000; HJ=false"
    "N=1000;  Nx=200;  Ny=200;  Npca=50;  Nz=10;  nR=10;  nP=1000; HJ=false"
    "N=1000;  Nx=200;  Ny=200;  Npca=50;  Nz=10;  nR=10;  nP=1000; HJ=true"
    "N=10000; Nx=2000; Ny=2000; Npca=200; Nz=0;   nR=1;   nP=1000; HJ=false"
    "N=10000; Nx=2000; Ny=2000; Npca=200; Nz=50;  nR=1;   nP=1000; HJ=false"
    "N=10000; Nx=2000; Ny=2000; Npca=200; Nz=50;  nR=1;   nP=1000; HJ=true"
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
  "try;${Param[i]};CCAinfMC;end;quit"

EOF

    echo qsub -N $ScriptNm $ScriptNm.sh
done

