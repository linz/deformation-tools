#!/bin/bash

rm out/*
echo "test 1 - read, stats" | tee -a out/tests.log
../gridtool read test1.grid  stats >> out/tests.log
echo "test 2 - zero" | tee -a out/tests.log
../gridtool test1.grid zero nearest_to test1_3rows.xy out/test1_zero.grid >> out/tests.log
echo "test 3 - smooth linear 1" | tee -a out/tests.log
../gridtool testlin1.grd smooth linear nearest_to testlin1.xy out/testlin1_out.grd >> out/tests.log
echo "test 4 - smooth linear 2" | tee -a out/tests.log
../gridtool test2.grid smooth linear nearest_to test1_3rows.xy out/test2_smoothl.grid >> out/tests.log
echo "test 5 - smooth quadratic" | tee -a out/tests.log
../gridtool test2.grid smooth quadratic nearest_to test1_3rows.xy out/test2_smoothq.grid >> out/tests.log
echo "test 6 - smooth cubic" | tee -a out/tests.log
../gridtool test2.grid smooth cubic nearest_to test1_3rows.xy out/test2_smoothc.grid >> out/tests.log
echo "test 7 - selection on wkt polygon 1" | tee -a out/tests.log
../gridtool test1.grid zero inside test.wkt out/test1_wkt_zero.grid >> out/tests.log
echo "test 8 - selection on wkt polygon 2" | tee -a out/tests.log
../gridtool test1.grid zero outside test.wkt out/test1_wkt_zero2.grid >> out/tests.log
echo "test 9 - multiply" | tee -a out/tests.log
../gridtool test1.grid multiply 3.0 out/test1_mult.grid >> out/tests.log
echo "test 10 - evaluate" | tee -a out/tests.log
../gridtool read test3.grid evaluate at test1.grid to out/test3_eval.txt >> out/tests.log
echo "test 11 - add" | tee -a out/tests.log
../gridtool test1.grid add test3.grid out/test1_add.grid >> out/tests.log
echo "test 12 - subtract" | tee -a out/tests.log
../gridtool test1.grid subtract test3.grid out/test1_subtract.grid >> out/tests.log
echo "test 13 - resize 1" | tee -a out/tests.log
../gridtool test1.grid resize 3 1 5 4 out/test1_resize1.grid >> out/tests.log
echo "test 14 - resize 2" | tee -a out/tests.log
../gridtool test1.grid resize -2 -3 10 11 out/test1_resize2.grid >> out/tests.log
echo "test 15 - trim" | tee -a out/tests.log
../gridtool test4.grid trim 1 out/test4_trim1.grid >> out/tests.log
echo "test 16 - extents_wkt" | tee -a out/tests.log
../gridtool read test1.grid extents_wkt out/gridextents.wkt >> out/tests.log
echo "test 17 - affected area" | tee -a out/tests.log
../gridtool read test1.grid affected_area inside affectedtest.wkt out/affectedout.wkt >> out/tests.log
echo "test 18 - selection with multiple where criteria 1" | tee -a out/tests.log
../gridtool test1.grid zero where v1 ">" 20 not v2 ">" 3.5 out/test1_where1.grid >> out/tests.log
echo "test 19 - selection with multiple where criteria 2" | tee -a out/tests.log
../gridtool test1.grid zero where v1 gt 20 and v2 le 3.5 out/test1_where2.grid >> out/tests.log
echo "test 20 - linzgrid with file based headers" | tee -a out/tests.log
../gridtool read test1.grid write_linzgrid NZGD2000 file header.txt out/test1_linzgrid.txt >> out/tests.log
echo "test 21 - linzgrid with inline headers" | tee -a out/tests.log
../gridtool read test1.grid write_linzgrid NZGD2000 "First line of header" "Second line of header" ThirdLineOfHeader resolution 0.001 out/test2_linzgrid.txt >> out/tests.log
echo "test 22 - linzgrid with columns" | tee -a out/tests.log
../gridtool read test1.grid write_linzgrid NZGD2000 "First line of header" "Second line of header" ThirdLineOfHeader columns v2+v1+v2 resolution 0.001 out/test3_linzgrid.txt >> out/tests.log
echo "test 23 - commands from file" | tee -a out/tests.log
../gridtool commands.txt >> out/tests.log
echo "test 24 - stats" | tee -a out/tests.log
../gridtool read csv test1.csv  stats >> out/tests.log
echo "test 25 - write no columns" | tee -a out/tests.log
../gridtool read csv test1.csv  write csv columns none out/test25.out >> out/tests.log
echo "test 26 - write selected columns" | tee -a out/tests.log
../gridtool read csv test1.csv  write csv columns v2 out/test26.out >> out/tests.log
echo "test 27 - write with selection" | tee -a out/tests.log
../gridtool read csv test1.csv  write out/test27.out where v1 gt 20 >> out/tests.log
echo "test 28 - add with selection" | tee -a out/tests.log
../gridtool test1.grid add test3.grid where v1 ">" 30 out/test28_add.grid >> out/tests.log
echo "test 28a - add with on_grid selection" | tee -a out/tests.log
../gridtool test1.grid add add1.grid where on_grid add1.grid out/test28a_add.grid >> out/tests.log
echo "test 28b - add with on_grid selection" | tee -a out/tests.log
../gridtool test1.grid add test3.grid where on_grid maxcols 1 add1.grid out/test28b_add.grid >> out/tests.log
echo "test 29 - replace" | tee -a out/tests.log
../gridtool test1.grid replace test3.grid out/test29_replace.grid >> out/tests.log
echo "test 30 - expand" | tee -a out/tests.log
../gridtool read test1.grid write columns none out/test30a.txt where nearest_to test30_pts.xy write columns none out/test30b.txt where nearest_to test30_pts.xy and expand 2.5 >> out/tests.log
echo "test 31 - align" | tee -a out/tests.log
../gridtool read align2.grid alignto align1.grid write out/test31.grid >> out/tests.log
../gridtool read out/test31.grid alignto align1.grid write out/test31a.grid >> out/tests.log
../gridtool read out/test31.grid resize relative -1 -1 1 1 alignto align1.grid write out/test31b.grid >> out/tests.log
../gridtool read out/test31.grid resize relative 1 1 -1 -1 alignto align1.grid write out/test31c.grid >> out/tests.log
echo "test 32 - trim" | tee -a out/tests.log
../gridtool read align2.grid resize relative -15 -15 45 25  trimto align1.grid write out/test32.grid >> out/tests.log
../gridtool read out/test32.grid resize relative -1 -1 1 1 trimto align1.grid write out/test32b.grid >> out/tests.log
../gridtool read out/test32.grid resize relative 1 1 -1 -1 trimto align1.grid write out/test32c.grid >> out/tests.log
../gridtool read align2.grid resize relative -15 -15 45 25  trimto buffer 2 align1.grid write out/test32d.grid >> out/tests.log
../gridtool read align2.grid resize relative -15 -15 45 25  trimto wkt test2.wkt write out/test32e.grid >> out/tests.log
../gridtool read align2.grid resize relative -15 -15 45 25  trimto buffer 3 wkt test2.wkt write out/test32f.grid >> out/tests.log

echo "test 33 - create" | tee -a out/tests.log
../gridtool create 25.0 29.0 0.25 -42.0 -41.0 0.1 columns de+dn write csv out/test33a.csv >> out/tests.log
../gridtool create extents align1.grid 0.25 0.1 columns de+dn write csv out/test33b.csv >> out/tests.log
../gridtool create extents wkt test.wkt 0.25 0.1 columns de+dn write csv out/test33c.csv >> out/tests.log

diff -q -r -B -b out check

