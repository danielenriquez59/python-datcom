# PyDATCOM worked examples

Runnable aircraft analyses that exercise the translated routines end to end.
Each script is self-contained and adds the repository root to `sys.path`, so
it runs from any working directory:

```bash
python examples/01_wing_tail_longitudinal.py
python examples/02_supersonic_wing_drag.py
python examples/03_datcom_input_file.py
```

Read `TRANSLATION_STATUS.md` for what is source-backed and `PHYSICS_REVIEW.md`
for the audit these examples were built to keep honest. **These are
demonstrations of the translated code paths, not validated engineering
predictions** — whole-aircraft results still depend on routines that remain
approximations.

## 01 — Conventional wing-tail longitudinal analysis

A 30 ft span straight-tapered wing with a conventional tail on a 25 ft arm.

Exercises: `WTGEOM` planform geometry for both surfaces, `INFTGM` downwash
synthesizing dimensions, the `DWASH` Section 4.4.1 downwash gradient, the
Section 4.3.1.2 lift carryover factors, `CLWBT` lift buildup, and `CMALPH`
zero-lift moment.

Produces an alpha sweep with downwash angle, tail local angle, and the wing
/ tail / total pitching moment breakdown, followed by `CL_alpha`,
`Cm_alpha`, `dCm/dCL`, static margin and neutral point.

This example is what caught two real defects: a per-radian tail lift slope
used against degree angles (a 57.3x error that produced a 3965% static
margin), and the aspect-ratio default described below. Broad sign-and-trend
unit tests had passed in both cases; only printing physical magnitudes
exposed them.

## 02 — Supersonic straight tapered wing drag

Exercises: `SUPDRG` zero-lift wave drag (Figure 4.1.5.2-58, sharp and round
leading edge), `SUPDRG` skin friction (exposed-MAC Reynolds number, source
roughness cutoff, Figure 4.1.5.1-27 Mach factor, Mach 3 cap), and `CORDSP`
supersonic section coordinates.

Produces a Mach sweep of wave drag, friction drag and their sum; a sharp
versus round leading-edge comparison; and a roughness sensitivity study.
The three CORDSP shapes each report a maximum full thickness of exactly the
requested `t/c`, which is the check for PHYSICS_REVIEW item 13.

## 03 — Original Digital DATCOM input deck

Parses `tests/fixtures/ex2.inp`, the AFFDL-TR-79-3032 Example Problem 2 deck,
and runs the translated pipeline across its flight-condition schedule.

Exercises: the namelist parser including the unindexed array continuations
used by the `MACH` and `ALSCHD` cards (PHYSICS_REVIEW item 15 — the deck's
eleven angles and four Mach numbers all parse), the state manager, and the
regime-dispatching `AerodynamicCalculator`.

Produces the parsed schedule, reference quantities, wing planform, regime
dispatch per Mach number, coefficients at each Mach, and a full angle sweep.

This deck is an exposed-wing-alone case, so it has no tail or body buildup.
Running it surfaced that the lift path fell back to a default aspect ratio of
6.0 when a deck supplies only raw `WGPLNF` geometry — this wing's exposed
aspect ratio is 1.80, so `CL` was roughly doubled. The lift path now resolves
the planform through `WTGEOM` first.

## Known wart

`pydatcom/aerodynamics/lift.py` treats `wing_aspect_ratio == 6.0` as a "not
supplied" sentinel, so a configuration that genuinely has an aspect ratio of
6.0 is silently recomputed from the planform. This is pre-existing behavior
and is left alone deliberately; it is recorded here and in
`tests/test_high_priority_physics.py` so it is not mistaken for intent.
