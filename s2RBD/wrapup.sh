for i in {1..60}
  do
#    cp -r model $i
#    sed -i 's/22028/'$RANDOM'/g' $i/mm_run*
#    cp  input/$i-s2.pdb  $i/s2-frbd-openmmawsem.pdb
    cd $i
#      sed -n '/Random/p' mm_run*
#      sbatch job.slurm
#      tail -1 trial-re.out
        tail -2107 relax.pdb >../re-lastf/end-$i.pdb
#      cp ../relax.slurm .
#      sbatch relax.slurm
#      tail -1 qc.data
#      cp ../ana.slurm .
#      sbatch ana.slurm
    cd ../


  done
