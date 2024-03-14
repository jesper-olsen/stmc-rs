marsaglia-rs
==============

Pseudo random number generator [1] based on Marsaglia's 1990 paper [2].
Produces floats in the interval [0;1) with a period of 2^144.

References
----------
[1] [The Art of Computer Programming, Donald E. Knuth, Vol. 2, Chapter 3](https://en.wikipedia.org/wiki/The_Art_of_Computer_Programming) <br/>
[2] "Toward a Universal Random Number Generator" by George Marsaglia, Arif Zaman and Wai Wan Tsang, Statistics & Probability Letters 8 (1990), pp. 35-39


Run
-----

```
% time cargo run --bin main --release

Frequency of samples in [0;0.4)
k, #samples, frequency, error
  1        2 0.50000000000000000000 +9.999999999999998e-2
  2        8 0.37500000000000000000 +2.5000000000000022e-2
  3       32 0.37500000000000000000 +2.5000000000000022e-2
  4      128 0.41406250000000000000 +1.4062499999999978e-2
  5      512 0.39843750000000000000 +1.5625000000000222e-3
  6     2048 0.39062500000000000000 +9.375000000000022e-3
  7     8192 0.39941406250000000000 +5.859375000000222e-4
  8    32768 0.39950561523437500000 +4.943847656250222e-4
  9   131072 0.40013885498046875000 +1.388549804687278e-4
 10   524288 0.40006065368652343750 +6.0653686523415296e-5
 11  2097152 0.39999675750732421875 +3.2424926758034545e-6

Histogram - 100 samples, 10 bins
Bin 0:    4
Bin 1:   18
Bin 2:   11
Bin 3:    8
Bin 4:   13
Bin 5:    5
Bin 6:    8
Bin 7:   11
Bin 8:    7
Bin 9:   15

Histogram - 100000 samples, 10 bins
Bin 0:        9998
Bin 1:       10012
Bin 2:       10005
Bin 3:        9994
Bin 4:        9998
Bin 5:        9995
Bin 6:       10000
Bin 7:       10002
Bin 8:        9996
Bin 9:       10000

min, max and their frequencies over 10000000000 samples
min 0; count 595
max 0.9999999403953552; count 597

real	0m57.333s
user	0m53.711s
sys	0m0.072s
```
