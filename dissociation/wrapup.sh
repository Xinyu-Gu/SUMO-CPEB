for i in {0..32}
  do
#    cp -r model $i
#    sed -i 's/E0=-44/E0=(-35+'$i')/g' $i/mm_run*
    cd $i
#      sed -n '/E0=/p' mm_run*
#      sbatch job.slurm
#      tail -1 trial.out
#      cp ../../CalcQ.slurm .
#      sbatch CalcQ.slurm
#      tail -1 amh.dat
#       cp ../ana.slurm .
#       sbatch ana.slurm
#       tail -1106 movie.pdb >../end-$i.pdb
#       tail -3313 movie.pdb >../end-$i.pdb
#      awk '{printf "%s\t%s\n", $1, $2}' interE.dat > dist.dat
      sed -n '202,$p' amh.dat >../$i-wham-2d.dat
#      cp d1_dist.dat ../d1/d1_$i.dat
##
    cd ../
#    awk '{printf "%s\t%s\n", $1, $2}' $i-wham-2d.dat > $i-wham.dat
#    echo $i-wham.dat $i 1 >>metadata
    echo $i-wham-2d.dat $i 0 1 0 >>metadata_2d
####
##
  done
###
#awk '{printf "%s %-1s %s\n", $1, $2-35, $3}' metadata > log
awk '{printf "%s %-1s %s %s %s\n", $1, $2-35, $3, $4, $5}' metadata_2d > log_2d
#mv log metadata
mv log_2d metadata_2d
