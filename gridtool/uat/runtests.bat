@echo off

setlocal
chdir /d %~dp0

del /q out\*.*
..\gridtool read test1.grid  stats >> out\tests.log
..\gridtool test1.grid zero test1_3rows.xy out/test1_zero.grid >> out\tests.log
..\gridtool testlin1.grd smooth linear testlin1.xy out/testlin1_out.grd >> out\tests.log
..\gridtool test2.grid smooth linear test1_3rows.xy out/test2_smoothl.grid >> out\tests.log
..\gridtool test2.grid smooth quadratic test1_3rows.xy out/test2_smoothq.grid >> out\tests.log
..\gridtool test2.grid smooth cubic test1_3rows.xy out/test2_smoothc.grid >> out\tests.log
..\gridtool test1.grid zero wkt:test.wkt out/test1_wkt_zero.grid >> out\tests.log
..\gridtool test1.grid zero outside wkt:test.wkt out/test1_wkt_zero2.grid >> out\tests.log
..\gridtool test1.grid multiply 3.0 out/test1_mult.grid >> out\tests.log
..\gridtool read test3.grid evaluate at test1.grid to out/test3_eval.txt >> out\tests.log
..\gridtool test1.grid add test3.grid out/test1_add.grid >> out\tests.log
..\gridtool test1.grid subtract test3.grid out/test1_subtract.grid >> out\tests.log
..\gridtool test1.grid resize 3 1 5 4 out/test1_resize1.grid >> out\tests.log
..\gridtool test1.grid resize -2 -3 10 11 out/test1_resize2.grid >> out\tests.log
..\gridtool test4.grid trim 1 out/test4_trim1.grid >> out\tests.log
..\gridtool read test1.grid extents_wkt out\gridextents.wkt >> out\tests.log
..\gridtool read test1.grid affected_area wkt:affectedtest.wkt out\affectedout.wkt >> out\tests.log
..\gridtool test1.grid zero where v1 ">" 20 and v2 "<=" 3.5 out/test1_where1.grid >> out\tests.log
..\gridtool test1.grid zero where v1 gt 20 or v2 le 3.5 out/test1_where2.grid >> out\tests.log
..\gridtool read test1.grid write_linzgrid NZGD2000 file header.txt out/test1_linzgrid.txt >> out\tests.log
..\gridtool read test1.grid write_linzgrid NZGD2000 "First line of header" "Second line of header" ThirdLineOfHeader resolution 0.001 out/test2_linzgrid.txt >> out\tests.log
..\gridtool commands.txt >> out\tests.log
