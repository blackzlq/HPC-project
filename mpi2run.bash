#!/bin/sh
#BSUB -J Liqin-project
#BSUB -o project_mpi2_output_file
#BSUB -e project_mpi2_error_file
#BSUB -n 21
#BSUB -R "span[ptile=7]"
#BSUB -q ht-10g
#BSUB cwd /home/zhang.liq/project/
work=/home/zhang.liq/project/
cd $work

tempfile1=hostlistrun
tempfile2=hostlist-tcp

echo $LSB_MCPU_HOST > $tempfile1
declare -a hosts
read -a hosts < ${tempfile1}
for ((i=0;i<${#hosts[@]};i+=2));
  do
   HOST=${hosts[$i]}
   CORE=${hosts[(($i+1))]}
	echo $HOST:$CORE >> $tempfile2
done

mpirun -np 21 -prot -TCP -lsf ./mpi2
