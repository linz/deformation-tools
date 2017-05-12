#!/bin/sh

echo "Okada test case 2" > out/run_okada_strike_test.out
../calc_okada -x -d -wl out/okada_2.wkt okada_2s.model point:2000:3000 out/okada_2s.out  >> out/run_okada_strike_test.out
python magnitude_angle.py out/okada_2s.out 0 2 9 >> out/run_okada_strike_test.out
../calc_okada -x -d -wl out/okada_2_30.wkt okada_2s_30.model point:232.051:3598.076 out/okada_2s_30.out >> out/run_okada_strike_test.out
python magnitude_angle.py out/okada_2s_30.out 0 2 9 >> out/run_okada_strike_test.out
