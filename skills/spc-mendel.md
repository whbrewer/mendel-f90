---
name: spc-mendel
description: Headless CLI workflow and parameter guidance for running the Mendel app in SPC, grounded in the Mendel’s Accountant literature and mendel.in semantics.
---

# SPC Mendel (Headless CLI)

Use this skill when running the Mendel app in SPC headless mode (no web UI), or when preparing biologically interpretable parameter sets for Mendel runs.

This skill provides:
- CLI workflows
- Parameter semantics
- Paper-grounded modeling intent
- Guardrails against biologically unrealistic configurations

It does **not** replace `help.html`; it explains how parameters behave in practice.

---

## Source of truth

Always defer to these artifacts in the following order:

1. Template defaults  
   - `apps/mendel/mendel.in`

2. Parameter definitions  
   - `apps/mendel/help.html`

3. Modeling literature (interpretive layer)  
   - ICCS 2007  
   - SCPE 2007  
   - BINP 2012 series  
   - TBMM 2015  
   - ICCS 2018  

---

## Headless workflow

### Submit a run

```bash
./spc submit mendel --params "pop_size=2000,num_generations=1000,mutn_rate=5.0"
```

### Start the scheduler

```bash
./spc scheduler
```

### Monitor output

```bash
./spc status <case_id>
cat user_data/cli/mendel/<case_id>/mendel.out
```

---

## REPL notes

```bash
spc> submit mendel pop_size=2000,num_generations=1000
spc> submit mendel --params "pop_size=2000,num_generations=1000"
```

---

## Parameter semantics (paper-grounded)

### mutn_rate
- New mutations per offspring
- Human-scale: ~3–10
- Viral-scale: ~0.1–1.0

### mutation_effect_distribution
- Natural distributions span orders of magnitude
- Uniform distributions are artificial

### fraction_beneficial
- Must be very small (<1e-3)

### selection_mode
- probability (default, realistic)
- truncation (artificial)

### environmental_variance / heritability
- Any noise weakens selection
- Heritability = 1.0 is artificial

---

## Model validity guardrails

The following combinations are biologically unrealistic:

- Uniform mutation effects
- Full truncation selection
- Zero environmental variance
- High beneficial fraction

---

## Paper-to-parameter intent mappings

### Mutation-count hypothesis (2012)

- Uniform mutation effects
- Truncation selection
- Near-zero noise

Expected: stabilization only under contrived conditions

### Selection threshold (STd)

Expected: linear accumulation below threshold

### Beneficial mutation threshold (STb)

Expected: beneficials do not rescue system

### TBMM (2015)

Expected: monotonic fitness decline

---

## Diagnostic expectations

- Fitness plateaus → artificial regime
- Stabilized mutation count → contrived selection
- Early beneficial dominance → expected

---

## Reproducing paper results

Record paper, figure, parameters, and seed for each run.

---

## Final notes

Mendel is a genetic accounting system, not an optimizer.
