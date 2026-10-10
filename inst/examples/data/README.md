# Synthetic E4 teaching data

`e4-baseline-1001.csv` and `e4-weak-1001.csv` are the prespecified first datasets
of the E4 main population study: N=160, four parallel items, group predictor,
and homogeneous independent Gaussian item noise with variance 1 (mean-scale T=.25).
They contain no OECD or other participant data. Both use P(group=a)=55/160.

Baseline: group difference 1.5 and stable within-group SD .7/.7.
Weak signal: total stable variance .01 and captured population variance .005.
`e4-provenance.csv` stores the generator seeds, E4 bootstrap seeds, population
targets and hashes of the source evidence. The generator and exporter are
`dev/e4_design.R` and `dev/prepare_e5_examples.R` in the source repository.

The installed `intervals_and_targets.R` lesson uses the same observations with
the ordinary default bootstrap seed 1. The E4 report used the recorded per-cell
seeds; its exact endpoints and failed-draw counts therefore need not equal the
lesson. Increasing B or changing the level leaves the original estimate intact.

The known generating quantities help interpret this example. They are not
inputs to `lambda_model()` and do not replace estimated measurement temperature.
