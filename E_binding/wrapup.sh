for i in {0..704}
 do
    mkdir sim-$i
    cd sim-$i
        mv ../sim-$i.fasta .
        cp ../6jxw-sim.pdb .
        python ~/CPEB3/analysis/Convert_openmm.py 6jxw-sim.pdb sim-$i.fasta B
        rm 6jxw-sim.pdb.bak
        mv 6jxw-sim.pdb sim-$i.pdb
        python ~/openawsem-master/mm_create_project.py sim-$i
        rm *clean*
        cp /home/xg23/CPEB3/SUMO/thread/model/* .

#        python mm_analysis.py sim-$i --thread 1 -t sim-$i-openmmawsem.pdb >trial.out

    cd ../

    mkdir complex-$i
    cd complex-$i
        cat ../sumo2.fasta ../sim-$i/sim-$i.fasta >complex-$i.fasta
        cp ../6jxw-openmmawsem.pdb .
        python ~/CPEB3/analysis/Convert_openmm.py 6jxw-openmmawsem.pdb complex-$i.fasta AB
        rm *.bak
        mv 6jxw-openmmawsem.pdb complex-$i.pdb
        python ~/openawsem-master/mm_create_project.py complex-$i
        rm *clean*
        cp /home/xg23/CPEB3/SUMO/thread/model/* .

#        python mm_analysis.py complex-$i --thread 1 -t complex-$i-openmmawsem.pdb >trial.out

    cd ../


 done
