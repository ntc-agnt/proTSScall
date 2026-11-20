# proTSScall — revision changelog (for the NTC ecosystem manager)

<!-- THIS FILE IS OWNED BY THIS REPO. The ecosystem manager (ntc-ecosystem-context) READS it to
track pipeline revisions and keep its overseer context current; it does NOT overwrite it.
Append a new entry at the TOP whenever you make a notable change (a revision / version bump /
results-affecting fix), so the manager has visibility. Keep entries to one or two lines. -->

Entry format (newest first):
```
## <YYYY-MM-DD> · <version/tag or claude/branch> · <NTC-### if any>
- <one-line summary>; affects results: <yes|no>
```

<!-- entries below, newest first -->

## 2026-11-20 · claude/vibrant-clarke-6uqmti · NTC-001
- add proTSScall_optimized.R (parallel to the original, which is untouched): Lever A vectorized dominant-TSS + dedup (byte-identical, proven by test/equivalence_dominant_dedup.R — 20k randomized trials, 0 mismatches, 2983 exercising the dedup tie-break) replacing the O(active×transcripts) per-row loops; NTC-001 guards (arg/input validation, exit-code checks on every external step, mkdir -p / ln -rsf, matrix-column assert); NTC-015 chr-naming assert before make_heatmap. affects results: no (byte-identical on valid input; new fail-fast on bad input/mismatch)
