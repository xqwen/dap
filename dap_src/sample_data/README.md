# Running sample data


All runs using 4 parallel thresds.

## Running with individual-level data


```
dap-g -d sim.1.sbamds.dat -t 4 -o output.sbams -l log.sbams
```

## Running with sufficient summary statistics

```
dap-g -d_est sim.1.est.dat -d_ld sim.1.LD.dat -d_n 343 -d_syy 515.6 -t 4 -o output.ss -l log.ss
```

## Running with z-scores

```
dap-g -d_z sim.1.zval.dat -d_ld sim.1.LD.dat -t 4 -o output.zval -l log.zval
```

