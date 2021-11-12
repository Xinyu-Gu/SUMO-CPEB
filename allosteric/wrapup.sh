for i in {0..70}
  do
#    cp -r model $i
#    sed -i 's/q0=0/q0=0.1+0.01*'$i'/g' $i/mm_run*
    cd $i
#      sed -n '/q0=/p' mm_run*
#      sbatch job.slurm
#      echo $i
      #tail -1 trial.out
       #cp ../ana.slurm .
       #sbatch ana.slurm
#      tail -1 info.dat
      sed -n '202,$p' info.dat >../$i-wham.dat
##
    cd ../
###
#    echo $i-wham.dat $i 200 >>metadata
###
#
  done
##
#awk '{printf "%s %-3s %s\n", $1, 1.15+$2*0.025, $3}' metadata > log
#mv log metadata
