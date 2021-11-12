for i in {0..70}
  do
    cd $i
      sed -n '202,$p' info.dat >../$i-wham.dat
#      cp d1_dist.dat ../d1/d1_$i.dat
##
    cd ../
    echo $i-wham-2dpc.dat $i 0 10000 0 >>metadata_pc
###
#
  done
##
awk '{printf "%s %-2s %s %s %s\n", $1, 0.1+0.01*$2, $3, $4, $5}' metadata_pc > log_2d
mv log_2d metadata_pc
