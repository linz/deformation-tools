#!/bin/sh

echo "Okada test case 2" > out/run_okada_test_case.out
../calc_okada -x -d -wl out/okdata_2.wkt okada_2s.model point:2000:3000 - >> out/run_okada_test_case.out
../calc_okada -x -d okada_2d.model point:2000:3000 - >> out/run_okada_test_case.out
../calc_okada -x -d okada_2t.model point:2000:3000 - >> out/run_okada_test_case.out


echo "-------- strike slip differentials" >> out/run_okada_test_case.out
../calc_okada -x -d okada_2s.model okada_test_point.dat out/okada_test_2s.out >> out/run_okada_test_case.out
python test_okada_differentials.py out/okada_test_2s.out >> out/run_okada_test_case.out
echo "-------- dip slip differentials" >> out/run_okada_test_case.out
../calc_okada -x -d okada_2d.model okada_test_point.dat out/okada_test_2d.out >> out/run_okada_test_case.out
python test_okada_differentials.py out/okada_test_2d.out >> out/run_okada_test_case.out
echo "-------- tensile differentials " >> out/run_okada_test_case.out
../calc_okada -x -d okada_2t.model okada_test_point.dat out/okada_test_2t.out >> out/run_okada_test_case.out
python test_okada_differentials.py out/okada_test_2t.out >> out/run_okada_test_case.out
