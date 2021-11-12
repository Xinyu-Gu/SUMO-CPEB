for i in {33..35}
  do
	for j in {0..8}
	  do
#	    cp -r model $i-$j
#	    cp mm_run.py $i-$j/
#	    sed -i 's/E0=-44/E0=(-35+'$i')/g' $i-$j/mm_run*
#            sed -i 's/r0=0/r0=(1+'$j'*0.5)/g' $i-$j/mm_run*
	    cd $i-$j
#	      sed -n '/E0=/p' mm_run*
#              sed -n '/r0=/p' mm_run*
#	      sbatch job.slurm
	      #tail -1 trial.out
	#      cp ../../CalcQ.slurm .
	#      sbatch CalcQ.slurm
#	      tail -1 amh.dat
#	       cp ../ana.slurm .
#	       sbatch ana.slurm
            sed -n '202,$p' amh.dat >../$i-$j-wham-2d.dat
	    cd ../
	#
    echo $i-$j-wham-2d.dat $i $j 1 10 >>log_2d
###
   done
#
  done
##
awk '{printf "%s %-1s %s %s %s\n", $1, $2-35, $3*0.5+1, $4, $5}' log_2d > log
cat metadata_2d log >log_2d
rm log	
mv log_2d metadata_2d
