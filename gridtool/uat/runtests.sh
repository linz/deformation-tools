#!/bin/bash

rm out/*
echo "test 1" >> out/tests.log
../gridtool read test1.grid  stats >> out/tests.log
echo "test 2" >> out/tests.log
../gridtool test1.grid zero nearest_to test1_3rows.xy out/test1_zero.grid >> out/tests.log
echo "test 3" >> out/tests.log
../gridtool testlin1.grd smooth linear nearest_to testlin1.xy out/testlin1_out.grd >> out/tests.log
echo "test 4" >> out/tests.log
../gridtool test2.grid smooth linear nearest_to test1_3rows.xy out/test2_smoothl.grid >> out/tests.log
echo "test 5" >> out/tests.log
../gridtool test2.grid smooth quadratic nearest_to test1_3rows.xy out/test2_smoothq.grid >> out/tests.log
echo "test 6" >> out/tests.log
../gridtool test2.grid smooth cubic nearest_to test1_3rows.xy out/test2_smoothc.grid >> out/tests.log
echo "test 7" >> out/tests.log
../gridtool test1.grid zero inside test.wkt out/test1_wkt_zero.grid >> out/tests.log
echo "test 8" >> out/tests.log
../gridtool test1.grid zero outside test.wkt out/test1_wkt_zero2.grid >> out/tests.log
echo "test 9" >> out/tests.log
../gridtool test1.grid multiply 3.0 out/test1_mult.grid >> out/tests.log
echo "test 10" >> out/tests.log
../gridtool read test3.grid evaluate at test1.grid to out/test3_eval.txt >> out/tests.log
echo "test 11" >> out/tests.log
../gridtool test1.grid add test3.grid out/test1_add.grid >> out/tests.log
echo "test 12" >> out/tests.log
../gridtool test1.grid subtract test3.grid out/test1_subtract.grid >> out/tests.log
echo "test 13" >> out/tests.log
../gridtool test1.grid resize 3 1 5 4 out/test1_resize1.grid >> out/tests.log
echo "test 14" >> out/tests.log
../gridtool test1.grid resize -2 -3 10 11 out/test1_resize2.grid >> out/tests.log
echo "test 15" >> out/tests.log
../gridtool test4.grid trim 1 out/test4_trim1.grid >> out/tests.log
echo "test 16" >> out/tests.log
../gridtool read test1.grid extents_wkt out/gridextents.wkt >> out/tests.log
echo "test 17" >> out/tests.log
../gridtool read test1.grid affected_area inside affectedtest.wkt out/affectedout.wkt >> out/tests.log
echo "test 18" >> out/tests.log
../gridtool test1.grid zero where v1 ">" 20 not v2 ">" 3.5 out/test1_where1.grid >> out/tests.log
echo "test 19" >> out/tests.log
../gridtool test1.grid zero where v1 gt 20 and v2 le 3.5 out/test1_where2.grid >> out/tests.log
echo "test 20" >> out/tests.log
../gridtool read test1.grid write_linzgrid NZGD2000 file header.txt out/test1_linzgrid.txt >> out/tests.log
echo "test 21" >> out/tests.log
../gridtool read test1.grid write_linzgrid NZGD2000 "First line of header" "Second line of header" ThirdLineOfHeader resolution 0.001 out/test2_linzgrid.txt >> out/tests.log
echo "test 22" >> out/tests.log
../gridtool read test1.grid write_linzgrid NZGD2000 "First line of header" "Second line of header" ThirdLineOfHeader columns v2+v1+v2 resolution 0.001 out/test3_linzgrid.txt >> out/tests.log
echo "test 23" >> out/tests.log
../gridtool commands.txt >> out/tests.log
echo "test 24" >> out/tests.log
../gridtool read csv test1.csv  stats >> out/tests.log
echo "test 25" >> out/tests.log
../gridtool read csv test1.csv  write csv columns none out/test25.out >> out/tests.log
echo "test 26" >> out/tests.log
../gridtool read csv test1.csv  write csv columns v2 out/test26.out >> out/tests.log
echo "test 27" >> out/tests.log
../gridtool read csv test1.csv  write out/test27.out where v1 gt 20 >> out/tests.log
echo "test 28" >> out/tests.log
../gridtool test1.grid add test3.grid where v1 ">" 30 out/test28_add.grid >> out/tests.log
