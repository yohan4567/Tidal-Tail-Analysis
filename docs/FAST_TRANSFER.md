# Cached-input transfers first

`prepare_fast_transfers.py` prioritizes NGC6397 and M4/NGC6121 independently of
the longer ten-cluster observation and calibration queues. Both must already
have complete observed selections and checksum-verified, uncut faint-star
calibration parents. The script does not accept a PM/CMD-selected science field
as a replacement calibration parent.

The small bright-star census is queried first at the Gaia AIP mirror, with ESA
as a fallback. The ADQL, complete-query checks and parent checksums are retained;
the responding endpoint is recorded. Changing the network endpoint does not
change the quality-support thresholds or scientific selection.

Each cluster independently goes through calibration, positive unit-mass
population checks, a small real-input response test, package export and a
separate relocated-path response test. Only the final
`conditional_model_package_ready` status identifies a deliverable run package.
These checks do not establish a validated total physical mass-loss rate.

The first preparation uses `outputs/fast_transfer_20261101_v1`:

- `model_preparation_status.json` records per-cluster readiness or failure.
- `packages/NGC6397_portable_v1.zip` and `packages/NGC6121_portable_v1.zip`
  are exported separately; an archive may exist before its final test finishes.
- `package_validation/<ID>/validation.json` records the required relocation test.
- Each archive contains `RUN_GUIDE.md`, pinned dependencies and input checksums.

NGC6397 passed that preparation. M4's initial bright queries failed transiently;
its completed retry census is verified separately and reused under
`outputs/fast_transfer_M4_20261101_v2`. That fresh snapshot uses gzip level 1
for intermediate CSV caches. This changes lossless file compression only;
the CSV values, model settings, source checks and validation steps are retained.
The M4 archive is under that version's `packages` directory, not the failed v1
attempt. Keep the recorded source and output versions separate.

Slow jobs continue separately. Do not run a second writer against the same
fast-transfer directory. Restart only after stopping the earlier writer and
reviewing any incomplete scientific outputs.
