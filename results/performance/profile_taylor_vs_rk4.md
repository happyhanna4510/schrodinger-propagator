# Taylor vs. RK4 profiling (dt=1e-5, tmax=0.01, K=4)

Test command:

```
./build/morse --evolve <method> --dt 1e-5 --tmax 0.01 --K 4 --profile
```

All runs used `N=2001`, `gamma=10`, `xmax=30`. Averages are reported per step over 1000 steps.

| Method | Variant | Step total [µs] | Series/work [µs] | RHS matvec [µs] | Reallocations |
|--------|---------|----------------:|-----------------:|----------------:|--------------:|
| Taylor | before  | 18.71 | 17.18 | 8.42 | 0.00 |
| Taylor | after   | 15.45 | 13.79 | 8.47 | 0.00 |
| RK4    | —       | 19.92 | 12.26 | 7.65 | 0.00 |

The optimized Taylor integrator now runs ~22% faster than RK4 (vs. ~6% slower before optimization) while keeping the same RHS evaluation count.
