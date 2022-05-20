# ClassGroups

:construction: Under construction :construction:

Python implementation of Class Groups of Imaginary Quadratic Fields, where elements are represented as positive definite binary quadratic forms.



## Benchmarks

Some timings of scalar multiplication using different algorithms. Nothng fancy, will eventually swap out with benchmarks

### Naive
```
Created group with 512 bits in 1.7881393432617188e-05 seconds
Performed a scalar multiplication with 511 bit secret in 0.1533510684967041 seconds
--------------------------------------------------
Created group with 1026 bits in 1.8835067749023438e-05 seconds
Performed a scalar multiplication with 1022 bit secret in 0.6584110260009766 seconds
--------------------------------------------------
Created group with 2049 bits in 3.695487976074219e-05 seconds
Performed a scalar multiplication with 2049 bit secret in 3.193856954574585 seconds
--------------------------------------------------
Created group with 512 bits in 1.2159347534179688e-05 seconds
Performed a scalar multiplication with 511 bit secret in 0.14811277389526367 seconds
--------------------------------------------------
Created group with 1024 bits in 2.002716064453125e-05 seconds
Performed a scalar multiplication with 1024 bit secret in 0.6391940116882324 seconds
--------------------------------------------------
Created group with 2047 bits in 3.814697265625e-05 seconds
Performed a scalar multiplication with 2046 bit secret in 3.183690071105957 seconds
--------------------------------------------------
```

### NUCOMP / NUDUPL
```
Created group with 512 bits in 1.71661376953125e-05 seconds
Performed a scalar multiplication with 511 bit secret in 0.12361025810241699 seconds
--------------------------------------------------
Created group with 1026 bits in 2.002716064453125e-05 seconds
Performed a scalar multiplication with 1026 bit secret in 0.474027156829834 seconds
--------------------------------------------------
Created group with 2049 bits in 4.100799560546875e-05 seconds
Performed a scalar multiplication with 2045 bit secret in 2.101064920425415 seconds
--------------------------------------------------
Created group with 512 bits in 1.0967254638671875e-05 seconds
Performed a scalar multiplication with 511 bit secret in 0.10851192474365234 seconds
--------------------------------------------------
Created group with 1024 bits in 1.9073486328125e-05 seconds
Performed a scalar multiplication with 1023 bit secret in 0.4572319984436035 seconds
--------------------------------------------------
Created group with 2047 bits in 3.981590270996094e-05 seconds
Performed a scalar multiplication with 2046 bit secret in 2.0988810062408447 seconds
--------------------------------------------------
```

### Gmpy2 + NUCOMP / NUDUPL

```
Created group with 512 bits in 3.790855407714844e-05 seconds
Performed a scalar multiplication with 507 bit secret in 0.0359799861907959 seconds
--------------------------------------------------
Created group with 1026 bits in 2.2172927856445312e-05 seconds
Performed a scalar multiplication with 1026 bit secret in 0.11339902877807617 seconds
--------------------------------------------------
Created group with 2049 bits in 3.409385681152344e-05 seconds
Performed a scalar multiplication with 2049 bit secret in 0.42668604850769043 seconds
--------------------------------------------------
Created group with 512 bits in 4.76837158203125e-06 seconds
Performed a scalar multiplication with 512 bit secret in 0.026320219039916992 seconds
--------------------------------------------------
Created group with 1024 bits in 5.0067901611328125e-06 seconds
Performed a scalar multiplication with 1021 bit secret in 0.10192584991455078 seconds
--------------------------------------------------
Created group with 2047 bits in 8.821487426757812e-06 seconds
Performed a scalar multiplication with 2046 bit secret in 0.41104888916015625 seconds
--------------------------------------------------
```