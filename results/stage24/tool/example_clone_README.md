# Runnable example — clone `L100615`

A real clone from the frozen WM989 benchmark, exported so the CLI can be run without any setup
beyond rebuilding the model artifact.

```text
clone_id     L100615
outer_fold   0
nuisance     log1p(n_naive_cells), log1p(n_naive1_cells), log1p(n_naive2_cells), log1p(n_naive3_cells)
             1.791759, 0.693147, 0.693147, 1.386294
```

Its frozen out-of-fold `pred_W5` scores, which the tool must reproduce when asked for
`--component fold0`:

```text
  Acid         0.6541259824971829   observed y = 0
  Cisplatin    0.24573579479469879   observed y = 0
  CoCl2        0.2127405206226653   observed y = 0
  Dabrafenib   0.24786048010260181   observed y = 0
  Doxorubicin  0.3928926208441807   observed y = 1
  Trametinib   0.1700245992727992   observed y = 0
```

Files: `example_clone_expression.npy` (36,601 CP10K/log1p values),
`example_clone_nuisance.txt` (the four-value nuisance block, comma separated).
