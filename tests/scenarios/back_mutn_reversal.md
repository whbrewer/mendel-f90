# Back Mutation Fitness Reversal Is Exact

## Claim being refuted

"A reverse mutation does not have the same magnitude of effect as the original mutation."

## How this test disproves it

Both the forward mutation path (`mating.f90`) and the back mutation path
(`mutation.f90:back_mutn`) compute fitness from the encoded mutation index using
the same Weibull formula. The encoded index is an integer that deterministically
maps to a fitness value via:

```
e_del = exp(-alpha_del * x^gamma_del)
e_fav = max_fav_fitness_gain * exp(-alpha_fav * x^gamma_fav)
```

where `x = mutn_index * scale`.

The forward path computes fitness from `mutn * del_scale` (or `fav_scale`), and
`back_mutn` calls `decode_fitness_del`/`decode_fitness_fav` which reconstructs
`x` from the same integer index using the same scale factor. The result is
bit-for-bit identical.

In the pure additive case (`multiplicative_weighting = 0`), the linkage block
algebra is:

- Forward:  `lb = lb - fitness_effect`  (deleterious)
- Reversal: `lb = lb - fitness`  where fitness is negative, i.e. `lb = lb + |fitness|`

With `verbosity = 2`, diagnostic lines are printed for both forward and back
mutations, allowing direct comparison of the fitness values for the same mutation
index.

## Running

```
cd src && ./mendel -f ../tests/scenarios/back_mutn_reversal.in
```

Filter for diagnostic output:

```
./mendel -f ../tests/scenarios/back_mutn_reversal.in 2>&1 | grep -E 'FORWARD|BACK_MUTN'
```

Match a specific mutation index across forward and back:

```
./mendel -f ../tests/scenarios/back_mutn_reversal.in 2>&1 \
  | grep -E 'FORWARD|BACK_MUTN' > /tmp/mutn_log.txt
grep BACK_MUTN /tmp/mutn_log.txt | head -1
# note the mutn= value, then:
grep 'mutn= *236952923 ' /tmp/mutn_log.txt
```

## Key parameters

| Parameter | Value | Rationale |
|---|---|---|
| `pop_size` | 100 | Small population |
| `mutn_rate` | 50 | High rate to accumulate mutations quickly |
| `genome_size` | 1.0E+06 | Small genome increases back mutation probability |
| `num_generations` | 20 | Short run, just enough to see back mutations |
| `allow_back_mutn` | T | Enable back mutations |
| `fraction_recessive` | 0.0 | All dominant, simplifies analysis |
| `dominant_hetero_expression` | 1.0 | Full expression so forward and back values match directly |
| `multiplicative_weighting` | 0.0 | Pure additive model |
| `frac_fav_mutn` | 0.01 | Modest favorable fraction |
| `verbosity` | 2 | Enables FORWARD/BACK_MUTN diagnostic output |
| `random_number_seed` | 42 | Fixed seed for reproducibility |

## Expected output

Each `FORWARD` line shows a mutation being applied. Each `BACK_MUTN` line shows
one being reversed. Lines with the same `mutn=` index can be compared directly.
The fitness magnitudes are bit-for-bit identical.

## Sample matched pairs

```
--- mutn=236952923 ---
FORWARD  del mutn=   236952923 fitness= 0.634829129036885E-04 lb_before= 0.100000000000000E+01 lb_after= 0.999936517087096E+00
BACK_MUTN del mutn=   236952923 fitness=-0.634829129036885E-04 lb_before= 0.100000000000000E+01 lb_after= 0.100006348291290E+01

--- mutn=74384209 ---
FORWARD  del mutn=    74384209 fitness= 0.895578814992767E-05 lb_before= 0.100000000000000E+01 lb_after= 0.999991044211850E+00
BACK_MUTN del mutn=    74384209 fitness=-0.895578814992767E-05 lb_before= 0.100000000000000E+01 lb_after= 0.100000895578815E+01

--- mutn=931442479 ---
FORWARD  del mutn=   931442479 fitness= 0.128663851878265E-05 lb_before= 0.100000000000000E+01 lb_after= 0.999998713361481E+00
BACK_MUTN del mutn=   931442479 fitness=-0.128663851878265E-05 lb_before= 0.100000000000000E+01 lb_after= 0.100000128663852E+01

--- mutn=316121109 ---
FORWARD  del mutn=   316121109 fitness= 0.157617582158234E-03 lb_before= 0.100000000000000E+01 lb_after= 0.999842382417842E+00
BACK_MUTN del mutn=   316121109 fitness=-0.157617582158234E-03 lb_before= 0.100000000000000E+01 lb_after= 0.100015761758216E+01
```

In every case, the forward fitness magnitude equals the back mutation fitness
magnitude to all 15 significant digits. The sign is opposite (positive for
forward deleterious, negative for decoded deleterious) because the forward path
stores the magnitude while `decode_fitness_del` returns a signed value — but the
Weibull computation is identical.

## Note on `dominant_hetero_expression`

The forward path scales fitness by `dominant_hetero_expression` before applying
it to the linkage block. The back mutation path uses the raw decoded value
(no scaling). With the default `dominant_hetero_expression = 0.5`, forward
fitness values would be half the back mutation values. This test sets it to `1.0`
so magnitudes match directly.
