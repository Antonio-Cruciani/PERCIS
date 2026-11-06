# PercIS with the source-only bag estimator

This separate C++ project is copied from
`/Users/antonio/Documents/GIT_CODE/PERCIS/code/PercIS/percis`.
The original folder and the copied legacy sources are unchanged.
`UPSTREAM_SHA256.json` records their original hashes.
The new executable is **`percis-source`**.

## Algorithm

The source-only version now follows PercIS's two-phase sample-size framework:

1. Estimate an upper bound on graph diameter using the original SCC/pivot
   method.
2. Draw a pilot sample. Estimate the average number of internal path vertices,
   the maximum mean score, and the maximum second moment.
3. Apply the PercIS Bennett sample-size formula to compute **one** production
   budget.
4. Discard the pilot and draw that many fresh production bags. Return only
   their mean scores.

There is no sequential error-based stopping, checkpoint schedule, or Hoeffding
fallback. The former `--max-samples` option and confidence-radius output have
been removed. Empirical Bernstein is still used inside the original PercIS
**pilot bound for average path length**; it is not a production stopping rule.

As in the supplied implementation, `-g N` bypasses the pilot and samples exactly
`N` production bags. Without `-g`, the pilot determines the fixed production
budget. Production and pilot counters are reported separately.

## Build

Requires C++17 and the header-only Boost.Multiprecision library. OpenMP is
optional. GCC and Boost are installed on this Mac:

```sh
cd /Users/antonio/Documents/codex/per_cent/tkdd/percis-source-only
make CXX=/opt/homebrew/bin/g++-16 OPENMP=1
make test CXX=/opt/homebrew/bin/g++-16 OPENMP=1
```

A serial build is available with `make`. Use `make clean` before changing
compiler or OpenMP settings. With OpenMP, the default uses all available
threads; `-t N` selects a specific number and `-t 0` selects the default.

## Run

Run the two-phase algorithm on the included Cit-HepPh data with eight threads:

```sh
./percis-source -d -t 8 \
  -o /tmp/source-centrality.csv \
  0.01 0.05 \
  examples/15_cit_hepph_rnd_init_50.txt examples/15_cit_hepph_el.txt
```

Here `epsilon=0.01` and `delta=0.05`. The default bag rule is the PercIS rule:
when there is more than one shortest path, sample
`max(1, floor(alpha * sigma_st))` paths, with `alpha=2.3`; otherwise sample one.
All inner draws are with replacement. `-s A` sets alpha explicitly and
`-s 0` (or a negative alpha) selects one path per pair. An unrepresentable bag
size raises an explicit error rather than overflowing.

Without an explicit `-s` or `--bag-size`, PercIS's pilot heuristic may switch
production to one path per pair when the pilot second-moment upper bound
exceeds `0.9` times the maximum-mean upper bound. The production variance bound
is adjusted for that change in estimator; pilot and production moments are
never mixed.

The fixed budget and fixed bag-size options remain available:

```sh
./percis-source -d -t 8 -g 100000 --bag-size 4 \
  -o /tmp/source-centrality.csv \
  --diagnostics /tmp/source-diagnostics.csv \
  0.01 0.05 \
  examples/15_cit_hepph_rnd_init_50.txt examples/15_cit_hepph_el.txt
```

This draws exactly 100,000 endpoint pairs, sampling four shortest paths for
each reachable pair. Unreachable pairs remain counted zero bags. `-g` counts
**bags**, never individual paths; `inner_paths` counts paths in reachable
production bags. The `--bag-size` override is mutually exclusive with `-s`.
With `-g`, epsilon and delta do not determine the sample count and there is no
accuracy certificate inferred from that budget.

`--pilot-samples N` overrides the initial pilot-size heuristic, with `N>=2`.
It cannot be combined with `-g`. The default heuristic is the supplied
PercIS expression `floor(max(1000, log(1/delta)/epsilon_normalized,
2*edge_diameter_bound*log(1/delta)))`.

A small runnable example with exact scores `(0,1/6,1/6,0)`:

```sh
./percis-source -d --indexed-ids -t 2 --bag-size 4 \
  0.01 0.05 examples/diamond.states examples/diamond.edges
```

Output is headerless `internal_id,estimate` CSV. Without `-o`, it goes to stdout;
progress and phase statistics go to stderr. `--diagnostics` writes
`internal_id,estimate,bag_sample_variance` with a header. This variance uses
production bag scores with denominator `m-1`; it is undefined (`nan`) for a
single production bag. The pilot does not contribute to either the estimated
centrality or its reported sample variance.

`--seed N` controls reproducibility. The same build, inputs, thread count,
production bag rule, production budget, and seed give the same production
output. The pilot uses a separate random stream. `-q` suppresses phase progress.
Use `--help` for the complete source-only CLI.

## Source-only estimator and range normalization

Let `S=sum_v x_v`. The estimated centrality is

```text
p_SO(v) = sum_{s != t, s != v, t != v} x_s * sigma_st(v)/sigma_st
          / ((n-2) * (S-x_v)).
```

Unreachable pairs contribute zero. Each bag draws its source with probability
`x_s/S`, then its target uniformly among the other vertices. For a bag of
`B_st` uniform shortest paths,

```text
d_v = (n-1)*S / ((n-2)*(S-x_v))
Y_v = d_v * (number of sampled paths containing v internally) / B_st.
```

Sources and targets are excluded. Both equal endpoint states and ascending
endpoint states are valid; no source–target state difference is used.
The factor `q(s,t)*d_v = x_s/((n-2)*(S-x_v))` gives unbiased bag scores.

To reuse PercIS's unit-range sample-size calculation, the pilot uses a
**single global normalization** `R=max_v d_v`: `Z_v=Y_v/R` and
`epsilon_normalized=epsilon/R`. First and second pilot moments are computed
from `Z_v`, while final estimates and bag variances are rescaled by `R` and
`R^2`. The path-length quantity rho remains an upper bound on `sum_v E[Z_v]`.

For normalized variance upper bound `v_hat`, maximum mean bound `p_hat`, and
path-length bound `rho_hat`, the production budget uses the original Bennett
objective

```text
sup_x log(2*rho_hat/(x*delta_production)) /
      (x*(1-x) * h(epsilon_normalized/(x*(1-x)))),
h(u) = (1+u)*log(1+u)-u,
delta_production = delta/2.
```

The domain is capped by the variance root, the maximum-mean bound, `1/2`, and
`rho_hat`. Pilot upper bounds reserve the remaining failure probability as in
the supplied code. This is the PercIS sample-size calculation, not a replacement
confidence-stopping method.

## Exact source-only reference in Julia

Use the exact source-only algorithm when evaluating `percis-source`; the
original source-target exact scores measure a different centrality.

```sh
julia --threads 8 exact_algorithm/run_source_only.jl -d \
  -o exact-source-only.csv \
  examples/15_cit_hepph_rnd_init_50.txt examples/15_cit_hepph_el.txt
```

The output has the same zero-based internal IDs and headerless CSV format as
the C++ results. `--plain` selects the original Julia one-score-per-line format.
The new exact program uses standard libraries only, with BigInt path counts
and BigFloat dependency accumulation. See [the exact algorithm instructions](exact_algorithm/README.md)
for the definition, normalization, and tests.

## Input and vertex IDs

- Edge files contain two nonnegative integer IDs per line. Duplicate edges and
  self-loops are ignored. The graph is undirected unless `-d` is supplied.
- State files contain one finite number in `[0,1]` per vertex. Blank lines and
  `#` comments, including trailing comments, are accepted.
- By default, IDs are assigned in **first-appearance order**, matching the
  supplied C++ reader and the included Cit-HepPh state-generation pipeline.
  A vertex appearing only in ignored self-loops is absent.
- `--indexed-ids` instead interprets state rows as IDs `0,...,n-1`; the state
  file defines the vertex set, including isolates. Do not use this option for
  the Cit-HepPh file with external paper IDs.
- `--id-map FILE` writes `internal_id,original_id` CSV.

All-node source-only centrality requires at least three vertices and at least
two positive states, ensuring `S-x_v>0` for every vertex. Equal positive states
are valid and yield normalized betweenness centrality. Undefined all-zero or
single-positive-state instances are rejected explicitly.

## Implementation details and corrections

The new code preserves the PercIS algorithmic framework and default bag policy;
it does not attempt to reproduce legacy bugs or produce bit-identical output:

- Production random streams do not replay pilot samples.
- Phase budgets are executed exactly, including with multiple threads.
- Zero-length bags enter the pilot length moments. An edge-diameter bound is
  converted to an internal-vertex bound by subtracting **one**, not two.
- The pilot length bound uses the sample-variance denominator `m-1`.
- A switch from bags to individual paths uses a valid single-path second-moment
  bound, rather than reusing the larger-bag second moment.
- The same Bennett objective is maximized over its full theorem domain using
  conservative interval upper bounds. The old epsilon-only cutoff and local
  grid scan could underestimate that supremum. This numerical correction can
  increase the computed budget relative to the supplied implementation.
- Bidirectional BFS expands whole frontiers from the side with smaller degree
  sum. Exact integer source weights, path counts, and path ranks prevent silent
  support loss or overflow when counts exceed `2^64`.
- Very small positive alpha values still produce at least one path per bag.

The source-only entry point is `source_main.cpp`, with implementation in
`src/SourceGraph.cpp`, `src/SourceSampler.cpp`, `src/SourceStructure.cpp`,
`src/SourceBounds.cpp`, and `src/SourceEstimate.cpp`. Tests cover exact scores,
bag variance, the two-phase sample accounting, bound calculations, diameter
bounds, parallel execution, and arbitrary-precision path sampling.

The original `main.cpp`, `src/Probabilistic.cpp`, and other legacy files remain
available for comparison. `make legacy CXX=/opt/homebrew/bin/g++-16` builds that
unchanged program. Legacy experimental CLI modes (top-k, VC/data-independent
ablations, uniform endpoint proposals, or comparison with external exact scores)
remain in that separate entry point; the source-only executable provides the
all-vertices PercIS run and fixed-budget run described above.
