for i in {1..20}
  do
    cp -r model $i
    sed -i 's/22028/'$RANDOM'/g' $i/mm_run*
    cp 1-frbd-openmmawsem.pdb $i/frbd-openmmawsem.pdb
 
    cd $i
      sed -n '/Random/p' mm_run*
    cd ../

  done


for i in {21..40}
  do
    cp -r model $i
    sed -i 's/22028/'$RANDOM'/g' $i/mm_run*
    cp 2-frbd-openmmawsem.pdb $i/frbd-openmmawsem.pdb

    cd $i
      sed -n '/Random/p' mm_run*
    cd ../

  done

for i in {41..60}
  do
    cp -r model $i
    sed -i 's/22028/'$RANDOM'/g' $i/mm_run*
    cp 3-frbd-openmmawsem.pdb $i/frbd-openmmawsem.pdb

    cd $i
      sed -n '/Random/p' mm_run*
    cd ../

  done

