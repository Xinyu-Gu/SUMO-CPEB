for i in {0..19}
do
mdconvert -o rrms-b-3spn2-$i.pdb -i $i rrms-b-3spn2.pdb
python mm_create_project.py rrms-b-3spn2-$i
mv clean.pdb clean-$i.pdb
done 

