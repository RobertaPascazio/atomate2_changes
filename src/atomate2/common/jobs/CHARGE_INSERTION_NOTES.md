# Known issue: `ChargeInterstitialGenerator` sometimes places candidate sites too close to existing atoms

Found 2026-10-01 through 2026-10-04 debugging ApproxNEB DFT reference runs for a
Zn-conductor screening effort (`zn_conductors/tier4`), while comparing DFT barriers
against CHGNet. Recorded here because it's the same code family as the
`charge_insertion_generator` threading added on this branch (`new-insertion-generator`),
and any real fix belongs here rather than as a one-off patch downstream.

## The bug

The charge-density-guided interstitial site generator (`ChargeInterstitialGenerator`,
from `pymatgen-analysis-defects`) occasionally proposes a candidate working-ion
insertion site that is unphysically close to an existing framework atom -- in the
worst observed cases, essentially coincident with it (as close as 0.058-0.129 A;
typical real bonding distances for the elements involved are 2.0-2.6 A). Checked
across 9+ independent endpoints in two materials (mp-1226810 CdGa2(SeS)2, mp-4809
Ga2HgS4) and found a bad contact in effectively all of them -- this is systematic for
these materials, not a rare edge case.

## Why it matters: two distinct failure modes downstream

1. **Loud failure.** The VASP ionic relaxation starting from the bad geometry goes
   numerically unstable (`eddrmm` -> `POTIM` -> `Positive energy` -> `brions`,
   repeating), and Custodian exhausts `max_errors_per_job` (default 5), FIZZLING the
   job. Annoying but at least visible.

2. **Silent corruption -- the more dangerous one.** The same instability can also
   happen to converge anyway, within VASP's SCF/force criteria, before Custodian's
   error budget runs out. The job then reports `COMPLETED` with no error at all, but
   the relaxation trajectory still dragged one or more *host framework* atoms to a
   wildly displaced, unphysical final position. Confirmed directly: mp-1226810's
   "endpoint 0" started with Zn only 0.853 A from a Ga atom, relaxed under ISIF=2
   (fixed cell -- see below), and reported success -- but a host Ga atom ended up
   displaced ~7.4 A from its starting position (out of an ~11-12 A cell), and a Cd
   atom by ~3.5 A. Nothing about the job's own state flags this; it only became
   visible downstream, when an ApproxNEB image interpolated between this corrupted
   endpoint and a good one produced a second host-atom collision (two Ga atoms on a
   collision course, worst at the image closest to the midpoint of the corrupted
   atom's bogus trajectory). The image-interpolation code itself was verified
   correct (does proper minimum-image/periodic-aware interpolation) -- the collision
   is a faithful reproduction of an already-bad input, not a second bug.

   Practical implication: **"COMPLETED" does not mean "physically valid" for an
   endpoint relaxation seeded by this generator.** Any analysis built on top of these
   results should not assume a successful VASP state implies a sane structure.

## Decision tree used to triage failures (MaxCorrectionsPerJobError only; a
`ValidationError` -- e.g. a corrupted/truncated `vasprun.xml` from a crash too fast
for Custodian's handlers to even engage -- is a separate failure class, not
necessarily reducible to a geometry problem the same way, and should be looked at
on its own)

- **Endpoint job FIZZLED**: check its own submitted structure directly for any
  unphysically short pairwise contact (not just involving the working ion). If
  found, that's the fix target -- move the offending atom before resubmitting.
- **Image job FIZZLED**: check both parent endpoints' own relaxations first, for
  whether either one displaced a *host* atom (excluding the working ion, which is
  expected to move a lot) by an implausible amount during its own relaxation. If
  yes, the image is just downstream of a corrupted endpoint -- nudging the image's
  geometry alone doesn't fix the real problem, and depending on how much you care
  about that particular hop, the pragmatic choice may be to drop it rather than
  chase a proper re-relaxation of the endpoint. If both endpoints look clean, the
  image's own interpolated geometry is the actual problem, and nudging the
  colliding atoms apart directly in that image is a reasonable fix.

A working implementation of this triage (`diagnose_fizzled_fw`) lives in
`Desktop/zn_conductors/tier4/mlip_dft_comparison.ipynb` as of 2026-10-04.

## What a real fix would need to do

Validate (or re-sample) candidate sites from `ChargeInterstitialGenerator` against a
minimum-distance criterion to every existing atom in the structure *before* handing
them off to relaxation -- i.e. catch this at generation time, not three VASP jobs and
a Custodian error cascade later. `atomate2`'s own ApproxNEB code already has a
`min_hop_distance` concept (`ApproxNebMaker(min_hop_distance=...)`, used to reject
hops that are too *short*); the missing piece is a symmetric check against sites that
are too *close* to an existing atom, applied where the candidate sites are first
proposed.

## Separate, related note: ISIF for ApproxNEB endpoint/image relaxations

Unrelated to the above, but found during the same debugging session: ApproxNEB
endpoint/image relaxations in this codebase default to `ISIF=3` (full cell relax --
ions, shape, and volume), inherited from `ApproxNebMaker` not overriding
`endpoint_relax_maker`/`image_relax_maker`. For comparison against force-field
calculators run with a fixed cell (`relax_cell=False`, the correct setting for
ApproxNEB-style calculations -- confirmed empirically that letting the cell relax
independently per endpoint breaks hop-matching/interpolation entirely, since
different endpoints then no longer share a common reference lattice), the DFT side
needs the matching `ISIF=2` (ions relax, cell fixed) applied specifically to the
endpoint/image jobs, e.g.:

```python
from atomate2.vasp.powerups import update_user_incar_settings
flow = update_user_incar_settings(flow, {"ISIF": 2}, name_filter="ApproxNEB image relax")
```

Host relaxation should stay at full `ISIF=3` (you want the true equilibrium host
lattice before building the supercell) -- only the per-endpoint/per-image steps need
pinning.
