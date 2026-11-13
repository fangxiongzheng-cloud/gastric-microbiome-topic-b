# Width-bin provenance (do not invent; do not claim false precision)

Recorded: 2026-11-13 04:00:05 America/Los_Angeles (UTC-8)
Box clock: 2026-11-13 20:00:05 Asia/Shanghai (UTC+8, box clock)

## Frozen addendum v0.6d (sha256 7de40171fcdbe274752591e79a820203a5e0f02142572d1b2b22baa44bb28303), lines 88–89

- Line width is split into **3 pre-fixed bins by |effect|**: `<0.2`, `0.2–0.4`, `>0.4`.
- v0.6d does **not** specify numeric point widths (pt) for those bins, and it does not specify a continuous-width formula.
- The later pictorial Figure 8 printed `0.55 + 1.80 × |effect|` pt. That continuous rule is a departure and must **not** be reused.

## Pre-result freeze-log note (2026-11-11 10:11 CST / 2026-11-10 18:11 America/Los_Angeles), line 68

- Repeats “width tiers |eff| <0.2, 0.2-0.4, >0.4”.
- Also does **not** state pt values.

## Original mechanical generator (scripts/s10_e11_schematic.py, preserved, not to be edited), line 39

```
def tier(x): x=abs(x); return 1.0 if x<0.2 else (2.5 if x<=0.4 else 4.5)
```

These matplotlib linewidths are **1.0, 2.5 and 4.5 pt**, with boundaries: |effect| < 0.2 → 1.0; 0.2 ≤ |effect| ≤ 0.4 → 2.5; |effect| > 0.4 → 4.5.

## How to use them

- Use the **v0.6d / freeze-log bin thresholds** as the frozen rule.
- Use the **original generator’s 1.0 / 2.5 / 4.5 pt values** as the restored implementation widths, citing that script line, not v0.6d, as the source of the point values.
- Do not invent a new pt scale. Do not claim that 1.0/2.5/4.5 pt appear in addendum v0.6d itself.
- Document this distinction in the figure QA and in any methods wording.

## E1a restored effect

- Restored display effect: sqrt(median of F2 and Korea stored noHP_rCLR R2) = 0.07157835860316486 (bin <0.2, 1.0 pt).
- Historical P<0.10 three-cohort effect (not for the new figure): 0.073948362138644086.
