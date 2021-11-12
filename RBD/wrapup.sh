for i in {1..60}
  do
    cd $i
#       sed -n '/Random/p' mm_run*
#	sbatch job.slurm
	tail -1 trial.out
    cd ../

  done

