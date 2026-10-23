# Migration Plan — coffea 0.7.31 → 2026.7.0 (BTVNanoCommissioning)

Branch: `update-coffea-2026` · Dev env: `btv_coffea_2026` (micromamba)

## 0. Strategy decision (the one that shapes everything)

**Target coffea 2026 in _eager_ / _virtual_ NanoEvents mode. Do NOT do a dask-awkward (delayed) rewrite.**

Rationale (from full codebase survey):
- All 18 workflow processors share an identical `process → process_shift` structure, **no shared base class**, plain-`dict` outputs, raw boolean-mask selections (no `PackedSelection`).
- The code is pervasively eager: `len(events)`, `np.ones(len(events))`, `.to_numpy()`, `events.caches[0]`, `np.full(len(weight), ...)` appear in every processor and in the histogrammer.
- Histograms are **already scikit-hep `hist.Hist`** (`import hist as Hist`); `coffea.hist` is only a comment reference in `plot_utils.py`. So the single biggest classic migration (coffea.hist → hist) is already done.

coffea 2026's `from_root` default mode is `'virtual'` (materialize-on-access) — arrays behave eager for `len()`/`.to_numpy()`. Running via `Runner` + non-dask executors keeps all the eager idioms working. A delayed rewrite would instead force changes in **every** processor + the histogrammer (dask-histogram, elimination of every `len()`/`to_numpy()`), for no physics benefit here. Eager confines the real work to ~2 files.

## 1. Blast radius

| Area | File(s) | Old (0.7) | Change |
|---|---|---|---|
| **Execution** | `runner.py` (only) | `run_uproot_job` + `iterative/futures/parsl/dask_executor` | ✅ Rewritten to `Runner` + executor classes (§3) |
| **JEC/JER build** | `utils/correction.py` (`JME_shifts`, ~1300–1669), `helpers/func.py` (factory) | `CorrectedJetsFactory.build(..., lazy_cache=events.caches[0])` | ✅ `lazy_cache` dropped; no pickle regeneration needed, that path is dormant (§4) |
| **Accumulators** | `helpers/func.py` (`dump_lumi`) | `column_accumulator`, `set_accumulator` | ✅ No change (both survive) |
| **Histogram fills** | `utils/histogramming/histogrammer.py` | eager flatten + numpy | ✅ Minimal; latent `if <ak> > <ak>` bug fixed (§6) |
| **Env/pins** | `setup.cfg`, `btv_env.yml` (was `test_env.yml`) | coffea 0.7.31 | DONE |
| **Outliers** | `BTA_producer.py`, `BTA_ttbar_producer.py`, `sf_ttdilep_kin.py` | extra `np.array(ak…)`, BDT `to_numpy` | ✅ Validated under eager (§8c, §8d) |

Everything else (corrections `extractor`/`LumiMask`/`Weights`/`BTagScaleFactor`, plotting scripts, condor submitters) works with minor/no changes — confirmed by introspection (§7).

## 2. Environment (done)
- `setup.cfg`: `coffea[dask,dask-awkward,xrootd]==2026.7.0`, `correctionlib>=2.6.0`, `python_requires=>=3.10`.
- `btv_env.yml` (renamed from `test_env.yml` in `6006cd1`; the CI `environment-file:` entries now follow the rename): python 3.10–3.12, coffea 2026.7.0 via pip, `parsl>=2024.12.09`, dropped `setuptools<=70.1.1` cap.
- Env `btv_coffea_2026` created (Python 3.12, numpy 2.x, awkward 2.11). Canonical CI env name stays `btv_coffea`.

## 3. Execution rewrite — `runner.py` (the core change)

`run_uproot_job(...)` and the `*_executor` **functions** are removed. Replacement is `coffea.processor.Runner` + executor **classes**.

Old:
```python
output = processor.run_uproot_job(
    sample_dict, treename="Events",
    processor_instance=processor_instance,
    executor=processor.futures_executor,
    executor_args={"skipbadfiles": ..., "schema": PFNanoAODSchema, "workers": N, "xrootdtimeout": 900},
    chunksize=args.chunk, maxchunks=args.max,
)
```
New:
```python
from coffea.processor import Runner, IterativeExecutor, FuturesExecutor, DaskExecutor, ParslExecutor
executor = FuturesExecutor(workers=args.workers)   # or IterativeExecutor() / DaskExecutor(client=client) / ParslExecutor()
runner = Runner(
    executor=executor,
    schema=PFNanoAODSchema,
    chunksize=args.chunk, maxchunks=args.max,
    skipbadfiles=args.skipbadfiles,
    xrootdtimeout=900,
)
output = runner(sample_dict, processor_instance, treename="Events")
```
Executor mapping (all confirmed present §7): `iterative`→`IterativeExecutor`, `futures`→`FuturesExecutor`, `dask/*`→`DaskExecutor(client=client)`, `parsl/*`→`ParslExecutor`. The parsl/dask cluster-config blocks stay; only the final `run_uproot_job` call site changes. These executors process chunks **eagerly** per worker (NOT the dask-awkward delayed graph — that's the separate `apply_to_fileset` path), which is why the eager-idiom processors keep working.

Watch-outs:
- **Fileset format**: `Runner.run` accepts old-style `{ds:[files]}` + `treename=`, so no fileset conversion needed for the `Runner` path.
- `skipbadfiles`/`xrootdtimeout`/`chunksize`/`maxchunks`/`schema` are now `Runner(...)` kwargs, not `executor_args`.
- Use non-dask executors (`IterativeExecutor`/`FuturesExecutor`) for first validation to guarantee materializing (virtual/eager) arrays.

## 4. JEC/JER — `correction.py` + `func.py`

`helpers/func.py:126-131` builds `JECStack`/`CorrectedJetsFactory` via `extractor().make_evaluator()`, pickled to a gzip blob.
`correction.py:1597-1601` applies:
```python
jets = correct_map["JME"]["jet_factory"][jecname].build(
    add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
    lazy_cache=events.caches[0],          # ← 0.7-only; events.caches gone; build() sig now build(injets)
)
met = correct_map["JME"]["met_factory"].build(events.PuppiMET, jets, {})
```
### ✅ DONE — and the blast radius is smaller than expected

Both call sites fixed (`correction.py:1597-1600`):
```python
jets = correct_map["JME"]["jet_factory"][jecname].build(
    add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll)   # lazy_cache dropped
)
met = correct_map["JME"]["met_factory"].build(events.PuppiMET, jets)  # 3rd arg dropped
```
Confirmed signatures in 2026.7.0: `CorrectedJetsFactory.build(self, injets)`, `CorrectedMETFactory.build(self, in_MET, in_corrected_jets)`.

**`helpers/func.py` needs NO changes** — all three constructors are unchanged:
`JECStack(corrections, jec=, junc=, jer=, jersf=)`, `CorrectedJetsFactory(name_map, jec_stack)`, `CorrectedMETFactory(name_map)`.

### ⚠️ KEY FINDING: the jetmet_tools/pickle path is DORMANT

Control flow in `JME_shifts`:
```
1042  if "JME" in correct_map:
1044      if "JME_cfg" in correct_map:   → correctionlib path  (LIVE)
1572      else:                           → jet_factory/met_factory pickle path (DEAD)
```
`JME_cfg` is set by both correctionlib branches; only the pickle branch
(`if "name" in conf["JME"].keys()`, `correction.py:370`) leaves it unset — and
**no campaign config in `AK4_parameters.py` contains a `"name"` key**. So:

- The **live JEC for every current campaign is correctionlib** (`jet_jerc.json.gz`, from `JME_path` for Run2-UL-NanoAODv15 or from `/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/<era>/latest/` for Run3). The migration to correctionlib JEC has **already happened**.
- The 6 tracked `jec_compiled*.pkl.gz` blobs are **stale artifacts that nothing loads**.
- **Pickle regeneration is NOT required.** (Verified they *would* fail: loading one under Python 3.12 raises `TypeError: code() argument 13 must be str, not int` — a cloudpickle code-object break. Rebuilding from raw `.txt` via `extractor`→`JECStack`→`CorrectedJetsFactory` works fine if ever needed.)
- The `.build()` fixes above are still correct and worth keeping — they de-fuse a latent coffea-0.7 call in the fallback path.

**Real remaining JEC work** is validating the *live* correctionlib block
(`correction.py:1044-1571`: `unc_jets`/`unc_met`, hand-rolled JES/JER sources, manual
`ak.flatten/unflatten/values_astype`) under awkward 2.11 / numpy 2.4 — this needs event data.
**✅ Resolved:** the block runs in every sweep (§8b, §8c) and in the CI commands with
`--isSyst all` (§8d).

### Validation result — corrections stack builds under coffea 2026

`load_SF(year, campaign)` run for every campaign in `AK4_parameters.py` (exercises
correctionlib, `lookup_tools.extractor`, `rochester_lookup`, `BTagScaleFactor`, `LumiMask`):

**11 / 13 OK** — all Run2 UL (`2016preVFP-UL`, `2016postVFP-UL`, `2017-UL`, `2018-UL`)
and all Run3 (`Summer22`, `Summer22EE`, `Summer23`, `Summer23BPix`, `Summer24`,
`Prompt25`, `prompt_dataMC`). Every successful campaign returns `JME_cfg` +
`JME_json_path`, **empirically confirming the correctionlib path is the live one** and no
campaign produces a `jet_factory` map.

**2 pre-existing failures (NOT coffea-related, unaffected by this migration):**
`Rereco17_94X` and `Winter22Run3` raise `KeyError: '<campaign>'` at
`correction.py:52` (`_cvmfs_dir` → `campaign_map()[campaign]`), reached from
`load_SF:106` during the **LUM** lookup — i.e. these legacy campaigns are simply absent
from `campaign_map()`. This aborts long before the JME branch, so the string-valued
`"JME"` entry for those two campaigns is latent/unreachable rather than the actual fault.

## 5. Accumulators — `helpers/func.py` `dump_lumi`
✅ **No change required** — `column_accumulator`, `set_accumulator`, and `accumulate` all survive in 2026.7.0 (§7). Optional cleanup: remove the dead `@property accumulator`/`self._accumulator` (never assigned) from all 18 processors when touching them. **✅ Done** (`31ef35e`, the 17 processors that had it).

## 6. Histogramming — `utils/histogramming/histogrammer.py`
- Fills already eager: `h.fill(temp_syst, flatten(genflavor), flatten(pruned_ev.SelJet[histname]), weight=flatten(...))`, `syst=np.full(len(weight), syst)`. Works under eager/virtual NanoEvents.
- `hist.Hist.fill` (hist 2.10) accepts the flattened awkward + numpy weights (API stable).
- **Latent bug to fix regardless**: `histogrammer.py:400-416` uses `pt if <ak_array> > <ak_array> else ...` — a Python `if` on a non-scalar awkward array. Rarely reached today; fix with `ak.where`. **✅ Fixed:** the `hl` record is built field by field with `ak.where` (and carries `charge`, §8c).

## 7. Ground-truth API check (installed coffea 2026.7.0, env `btv_coffea_2026`)

Stack: coffea **2026.7.0**, awkward **2.11.0**, hist **2.10.1**, uproot **5.7.5**, correctionlib **2.9.0**, numpy **2.4.6**, dask-awkward **2026.2.1**, Python **3.12.13**.

**GONE (all in `runner.py` only):** `run_uproot_job`, `iterative_executor`, `futures_executor`, `parsl_executor`, `dask_executor`.

**Replacements ALL present:** `Runner`, `IterativeExecutor`, `FuturesExecutor`, `DaskExecutor`, `ParslExecutor` (only `WorkQueueExecutor` absent — unused here).

**Accumulators ALL present (no change):** `accumulate`, `column_accumulator`, `set_accumulator`, `dict_accumulator`, `defaultdict_accumulator`, `list_accumulator`, `ProcessorABC`. Schemas `PFNanoAODSchema`/`NanoAODSchema`/`NanoEventsFactory` present.

**Corrections/lookup ALL present:** `extractor`, `txt_converters`, `rochester_lookup`, `dense_lookup`, `LumiMask`, `Weights` (`Weights(size, storeIndividual=False)`), `PackedSelection`, `JECStack`, `CorrectedJetsFactory`, `CorrectedMETFactory`, `BTagScaleFactor`. Alt modern path also present: `dataset_tools.apply_to_fileset/preprocess/max_chunks/slice_chunks`.

**CHANGED signatures (must edit call sites):**
- `CorrectedJetsFactory.build(self, injets)` — `lazy_cache` removed; `events.caches` gone.
- `Runner.__init__(executor, pre_executor=None, chunksize=100000, maxchunks=None, skipbadfiles=False, xrootdtimeout=60, schema=NanoAODSchema, format='root', ...)`; `Runner.run(fileset, processor_instance, *, treename=None, ...)`.
- `NanoEventsFactory.from_root(file, *, mode='virtual', treepath=…, schemaclass=NanoAODSchema, …)` — default mode now **'virtual'**.

**Net effect:** only **two call-site signatures** actually change (`run_uproot_job→Runner` in `runner.py`; drop `lazy_cache` in `correction.py`). Low-friction migration, not a rewrite.

## 8. Peripheral / infra
- Execution API is **only** in `runner.py` (verified) — no other file uses `run_uproot_job`.
- Env naming: canonical name stays `btv_coffea` (CI `activate-environment: btv_coffea`; hardcoded in runner.py:722, condor_lxplus/submitter.py:231/235, scripts/fetch.py). Keep the migrated production env named `btv_coffea`; `btv_coffea_2026` is a dev sandbox only.
- Postprocessing scripts (`scripts/plotdataMC.py`, `make_template.py`, `comparison.py`, `plotSysts.py`, …) load `.coffea` via `coffea.util.load`; output object shape unchanged in eager mode → expected to keep working (watch `hist` version drift in `plot_utils.py`).
- CI workflows (`.github/workflows/*_workflow.yml`, `validation.yml`) build `btv_env.yml` and run one workflow each — these become the end-to-end migration gate. CI pins Python 3.10, while the sweeps ran on 3.12.

## 8b. END-TO-END VALIDATION ✅ (P1+P2+P3 verified with real data)

```
python runner.py --wf QG_dijet --json metadata/Summer24/MC_Summer24_2024_QG_dijet.json \
       --campaign Summer24 --year 2024 --overwrite --max 1 --limit 1
```
**Result: PASSES under coffea 2026.7.0** (executor `futures` / `FuturesExecutor`, 3 workers)
```
Dataset validation complete. 11 valid samples remaining.
Preprocessing 100% 11/11 [0:00:20]
Processing    100% 11/11 [0:02:24]
Saving output to hists_QG_dijet_MC_Summer24_2024_QG_dijet.coffea   (97 kB)
```
Output verified by reloading with `coffea.util.load`: **11 datasets × 45 `hist.Hist`**,
plus `fname`/`run`/`lumi`/`sumw` (confirms `column_accumulator`/`set_accumulator` survive).
`sumw = 7.44e11`; filled histograms carry sensible weighted sums (~1.07e11 each).

> Note: only the `ObjSelJet_*` histograms are filled (9/45 per dataset). This is **not** a
> migration issue — in the current working tree `dijet.py` has `CenJet`/`FwdJet`/`RndJet`/
> `LeadJet`/`SubleadJet` commented out (uncommitted local edits), so only `SelJet` is
> defined. The `qgtag` collection books histograms for all object categories regardless.
> **✅ Resolved:** the categories are back in `dijet.py`; the Oct 2026 re-sweep fills 396/495.

### Runtime fixes required (7 iterations)

| # | Site | Cause | Category |
|---|---|---|---|
| 1 | 37 sites / 10 files | `events.Jet = …` attribute field assignment | **awkward 2** |
| 2 | `correction.py` ×4 | field named `rho` clashes with vector's azimuthal-radial alias vs `pt` → renamed `event_rho` | **vector 1.8** |
| 3 | `get_corr_inputs:645` | correctionlib input-name map `Rho`→`rho` had to follow the rename → `event_rho` | knock-on |
| 4 | `MuonScaRe.py:178` | `random.seed(np.uint32)` — Py3.11 removed the `hash()` fallback → `int(seed)` | **Python 3.11+** |
| 5 | `MuonScaRe.py:242` | masked in-place assign `k_f[cond] = …` → `np.where` | **awkward 2** |
| 6 | `correction.py:1086` | `GenJet[0]` on events with **zero** GenJets → pad to length ≥1 | **latent bug** |

### Lessons vs. the original risk ranking (§10)
- The **core strategy bet was correct**: eager mode via `Runner` worked from the first run and
  never regressed. Every failure was in the processor/corrections layer (the predicted P3 class).
- **numpy 2.x — ranked risk #1 — caused zero failures.**
- Real cost was **awkward-2 strictness** (3 of 6 fixes), plus two categories the plan never listed:
  **Python-version drift** (coffea 2026 forces ≥3.10; the solve picked 3.12) and a **latent data bug**.
- Fix #6 matters beyond this migration: with a non-zero index the old code would have silently
  produced *wrong* `Genpt` values instead of crashing.

## 8c. P3 WORKFLOW SWEEP RESULTS (41 / 41 validated end-to-end)

All runs: `--max 1 --limit 1` on Summer24 MC unless noted. Env `btv_coffea_2026`
(coffea 2026.7.0). The registry grew from 33 to 41 workflows with the
`QG_photondijet*` / `QG_trijet*` additions. Those 8 variants (`cebc2b3`) now live on
the `split/qg-sf` branch; `update-coffea-2026` registers 33.

**PASS (36).** Sweeps 1-2 (23): `QG_dijet` `QG_DY` `QG_photonjet` `QG_zerobias`
`QG_pfjet` `ctag_Wc_sf` `ectag_Wc_sf` `ctag_Wc_noMuVeto_sf` `ctag_Wc_WP_sf`
`ectag_Wc_WP_sf` `ctag_DY_sf` `ectag_DY_sf` `DY_sfl` `eDY_sfl` `ctag_ttsemilep_sf`
`ectag_ttsemilep_sf` `ctag_ttsemilep_noMuVeto_sf` `ttsemilep_sf` `c_ttsemilep_sf`
`sf_ttsemilep_tnp` `QCD_sf` `example` `validation`.
Sweep 3 (13): `QG_photondijet` `QG_photondijet_quark` `QG_photondijet_softprobe`
`QG_photondijet_softprobe_quark` `QG_trijet` `QG_trijet_gluon` `QG_trijet_dijetave`
`QG_trijet_dijetave_gluon` `ttdilep_sf` `ctag_ttdilep_sf` `ectag_ttdilep_sf`
`emctag_ttdilep_sf` `sf_ttdilep_kin`.

Sweep 3 follow-up (5): `BTA` `BTA_addPFMuons` `BTA_addAllTracks` `BTA_ttbar`
`QCD_smu_sf` — these five failed on first contact and pass after the fixes in
"Fixing the sweep-3 failures" below. **No workflow in the registry is now unvalidated.**

Every sweep-3 PASS was checked for **non-vacuity**: the `.coffea` was reloaded with
`coffea.util.load` and at least one `hist.Hist` verified to have a non-zero weighted sum
(e.g. `emctag_ttdilep_sf` 36/58 histograms filled, `sumw` 4.05e6; the four
`QG_photondijet*` 540/615 filled each). An exit code of 0 is not sufficient evidence —
see the BTA finding.

### Sweep 3 (Aug 2026): the 8 new QG variants + the 10 previously-blocked workflows

Missing filesets were covered by reusing sibling JSONs rather than fetching from DAS:
`QG_photondijet*` → `MC_Summer24_2024_QG_photonjet.json`, `QG_trijet*` →
`MC_Summer24_2024_QG_dijet.json`, all ttdilep variants and `BTA_ttbar` →
`MC_Summer24_2024_ctag_ttsemilep_sf.json` with `--only TTto2L2Nu…`.

**The BTA skip-guard was bypassed without editing source**, by running with
`--isSyst JP_MC`: `BTA_producer.py:66-72` prefixes `systematic/` to the `gfal-ls` path
when `self.isSyst` is truthy, and that path does not exist centrally, so the early
`return` does not fire.

> **§8c's suspicion is now confirmed: the BTA CI gate was passing vacuously.** With the
> guard bypassed, all three `BTA*` producers crash immediately under coffea 2026. The gate
> had been green only because the processor returned before doing any work.

### Migration bugs found by sweep 3
| Fix | File | Cause |
|---|---|---|
| `PtEtaPhiECandidate` missing `charge` | `utils/histogramming/histogrammer.py:419` | the `hl` ("harder lepton") record is zipped from pt/eta/phi/energy only; coffea 2026 validates a behaviour's required fields **at construction**, where 0.7 did not. Broke every workflow booking `hl`, i.e. `emctag_ttdilep_sf`. **Fixed** by zipping `charge` with the same `ak.where(_mu_is_harder, …)` as the other four fields; `emctag_ttdilep_sf` then passes with `hl_ptratio` filled. |

That is the **only** migration bug in 18 runs — consistent with sweep 2 (13 workflows,
0 bugs). The migration itself is converged.
*(Superseded by §8d: the CI test file reached three more awkward-2 bugs in paths these
Summer24 inputs never exercise. All three are fixed.)*

### Fixing the sweep-3 failures

The three reported defects were each **masking further ones**: they aborted the
processor before code that had never run under coffea 2026 could be reached. Fixing
one exposed the next, so the five workflows took six fixes in total. Every fix below
was verified by re-running the workflow on real data.

| # | Site | Cause | Category |
|---|---|---|---|
| 1 | `BTA_producer.py` (genpt) | `events.GenJet[genJetIdx]` with the invalid-match index clamped to `0`, which is *itself* out of range in events with zero GenJets. Padded with `ak.pad_none(..., 1)`. Second instance of §8b fix #6. | **awkward 2** |
| 2 | `BTA_producer.py`, `BTA_ttbar_producer.py` (jet ID columns) | `jet.jetId` unguarded, but **NanoAODv13+/v15 dropped `Jet_jetId`** (confirmed absent from the Summer24 v15 branch list). `utils/selection.py` already guarded the same access; the BTA producers never got it. | **NanoAOD version** |
| 3 | `QCD_soft_mu_validation.py:179` | `Jet[:, 1]` while the selection requires only `>=1` jet. Padded so 1-jet events yield `-1`; acceptance unchanged. | **latent bug** |
| 4 | `BTA_producer.py` ×4 (`sign2D`, `sign3D`, `ip2/3dsign_jetref`) | `ThreeVector.dot(LorentzVector)` — vector>=1.8 refuses mismatched dimensionality instead of silently projecting. Take `jet.pvec`. | **vector 1.8** |
| 5 | `BTA_producer.py` ×2 (`ptrel`) | `cross` on two 4-vectors — "cross is only defined for 3D vectors". Take the spatial parts. | **vector 1.8** |
| 6 | `BTA_producer.py:979` | a field literally named **`deltaR`** on a Lorentz-vector record. vector>=1.x defines `deltaR` as a *method*, which shadows the field, so the later read returned a bound method (`unknown type method`). Renamed the intermediate field to `deltaR_jet`; the output branch is still `deltaR`. | **vector 1.8** |
| 7 | `QCD_soft_mu_validation.py` (JetSVs) | `ak.pad_none(JetSVs.jetIdx, len(filtered_events.Jet))` — `len` of a jagged array is the *event* count, not the per-event jet count, so the pad target was meaningless. Removed: the `valid_indices` filter immediately below already does what the pad was reaching for. | **latent bug** |

Fix 2 needed a small, behaviour-preserving refactor: `utils/selection.py` now exposes
`jet_id_mask(events, campaign)` (the per-jet ID with no kinematic cuts) and `jet_id`
is that mask plus the pt/eta cuts, exactly as before. `helpers/BTA_helper.py` gained
`jet_id_columns(...)`, shared by both BTA producers. Where `Jet_jetId` is absent the
framework only implements TightLepVeto, so all three ID columns are filled with it —
TightLepVeto implies Tight implies Loose, so a flagged jet genuinely passes that
column, and jets passing only a looser ID are conservatively `False` rather than
guessed. **This is a physics-content change to the BTA ntuple on v15 inputs and wants
an owner's eye.**

> **§8c's suspicion is confirmed twice over.** The BTA gate was not merely vacuous —
> the code behind it had *six* separate coffea-2026 incompatibilities, none of which
> any green CI run could have caught. `grep -c "skip "` on the fixed runs is 0, and
> the logs show the producers reaching their ntuple write (which then fails only on
> central-EOS permissions), i.e. full execution.

Fix 6 generalises: **any field whose name collides with a `vector` method now
resolves to the method.** This is the same failure class as the `rho` → `event_rho`
rename (§8b fix #2) and is worth a repo-wide grep before the next campaign.
**✅ Done** (§9 P5): no live clash.

### Caveats on this table
- *(Resolved Oct 2026, see §9 P4: all 23 re-run and passing.)* The 23 workflows from sweeps 1-2 were **not** re-run in sweep 3, so they are not verified
  against the later `PackedSelection` refactors (`8983295`, `9733c9c`, `240ba5f`), the new
  awkward interfaces / explicit TTree writing (`6fa8cda`), or the array-writer output change
  (`0e4b8ed`).
- `QCD_smu_sf` and the four `BTA*` workflows were re-verified after the fixes, together
  with `example`, `ttsemilep_sf` and `QG_trijet` as regression checks on the
  `jet_id` / `jet_id_mask` refactor (identical filled-histogram counts to sweep 3).
- Sweep 3's first pass raced a concurrent edit adding `ttbar_reweights` to every processor
  `__init__`; 6 workflows aborted with
  `TypeError: NanoProcessor.__init__() got an unexpected keyword argument 'ttbar_reweights'`.
  Those runs were repeated against the settled tree and are **not** counted as failures. If
  that keyword is ever added to `runner.py` again ahead of the processors, every workflow
  fails at construction — worth a CI smoke test that merely instantiates all 41 processors.
  **✅ Done:** `scripts/check_workflows.py` + `check_workflows.yml` (§9 P4).

## 8d. CI COMMANDS ON PYTHON 3.10 (Oct 2026)

Every distinct `runner.py` command in `.github/workflows/*.yml` (20) was run locally in a
fresh env built from `btv_env.yml` with `python=3.10` (3.10.22, coffea 2026.7.0), on
`metadata/test_bta_run3.json` (Summer23 PFNano), with the same non-vacuity check as §8c.

**All 20 pass after three fixes.** None of the failures were Python-3.10 specific. All three
were coffea-2026 / awkward-2 bugs in code paths that the Summer24 BTV Nano used by the sweeps
never reaches, because those inputs have no `PFCands`, no `JetSVs` and never zero-GenJet MC
events in the BTA ttbar sample. The CI test file has all three, so CI would have caught them on
the PR. Upstream master's CI passes the same commands under coffea 0.7.

| # | Site | Cause | Hit by |
|---|---|---|---|
| 1 | `BTA_ttbar_producer.py` (genpt) | §8c BTA fix 1 never ported: placeholder GenJet index 0 is out of range in events with zero GenJets. Same `ak.pad_none(events.GenJet, 1)` fix. | `BTA_ttbar` |
| 2 | `helpers/func.py` `PFCand_link` | multi-jet branch returned a `dict` with **int** keys, assigned onto `pruned_ev`. awkward 2 rejects it (`keys of dicts in 'fromiter' must all be strings`). Now string keys, `ak.zip(..., depth_limit=1)`. Nothing downstream reads the per-jet form (the histogrammer needs `PFCands.pt`, the array writer drops `PFCands`). | `ttsemilep_sf`, `c_ttsemilep_sf`, `ttdilep_sf`, `validation` |
| 3 | `QCD_validation.py` (JetSVs) | §8c fix 7 never ported from the soft-mu sibling: attribute field assignment plus a pad sized by `len()` (the event count). | `QCD_sf` |

The three `BTA*` CI jobs remain vacuous for `BTA`/`BTA_addPFMuons`/`BTA_addAllTracks`: their
outputs already exist on central EOS, so the producer skips every file. `BTA_ttbar` has no
central output yet and runs fully.

Not covered locally: lxplus mounts all of `/cvmfs`, but the CI workflows mount only
`grid.cern.ch`, while the corrections read from `/cvmfs/cms-griddata.cern.ch`. Master has the
same mount list and passes, so this is probably fine, but it is the one CI difference this run
could not reproduce.

## 9. Phasing / PR breakdown
- **P0 (done)**: pins + env + `import coffea` smoke test.
- **P1 (done)**: `runner.py` execution rewrite to `Runner` + executor classes (§3, §8b).
- **P2 (done)**: JEC/JER — `lazy_cache` / `events.caches[0]` dropped, `.build()` calls fixed. Pickle regeneration is not needed: that path is dormant (§4). The live correctionlib `JME_shifts` block runs in every sweep.
- **P3 (done)**: all 41 registered workflows validated end-to-end (§8c). This branch now registers 33; the 8 `QG_photondijet*`/`QG_trijet*` variants live on `split/qg-sf`. The numpy-2 alias audit is clean: no `np.float_`/`np.bool8`/`np.NaN`/`np.int0`/`np.complex_`/`np.unicode_`/`np.string_` left in `src/`, `scripts/` or `runner.py`.
- **P4 (partly done)**:
  - ✅ parsl/dask cluster configs build under parsl 2026 / dask-jobqueue 0.9 (`43eac31`). Real SLURM/HTCondor submission through them is **not** verified.
  - ✅ `condor/submitter.py` works from an EOS checkout (`96b2bab`, `561500f`) and is used for production on lxplus.
  - ✅ CI `environment-file:` entries point at `btv_env.yml` again. They still named `test_env.yml` after the rename, so every CI job failed at environment setup.
  - ✅ All 20 CI `runner.py` commands pass locally on Python 3.10 after three fixes (§8d). Still to see: the real GitHub run on the PR to master. CI only triggers on `master` pushes and PRs.
  - ✅ Re-ran the 23 sweep-1/2 workflows (Oct 2026, same settings as §8c, iterative executor) against the later PackedSelection / TTree-writing / array-writer changes and `3c2b0ec`, `7f2e3a1`: **23 / 23 pass with filled histograms**. `sf_ttsemilep_tnp` needs `--selectionModifier tt_semilep_mu|tt_semilep_el` since `f6f4bd4`, and passes in both channels (603/800 and 355/800 filled).
  - ✅ CI smoke test: `scripts/check_workflows.py` builds and pickles every registered workflow through `runner.build_processor` (the code path a run uses), with no events, proxy or CVMFS. It runs in `.github/workflows/check_workflows.yml`. `sf_ttsemilep_tnp` is built once per channel because it rejects an empty `--selectionModifier` (since `f6f4bd4`).
- **P5 (started)**:
  - ✅ `plot_utils.py` rebin helper uses `h.storage_type()` instead of `_storage_type()`, which hist 2.10 deprecates with a `FutureWarning`.
  - ⬜ Run the post-processing scripts (`plotdataMC.py`, `make_template.py`, `comparison.py`, …) on a coffea-2026 `.coffea` output.
  - ✅ Removed the unused `accumulator` property from the 17 processors that carried it (§5). `ProcessorABC` only requires `process`.
  - ✅ Checked for record fields whose names collide with `vector` methods (§8c fix 6): every field name assigned in `src/` was compared with the attributes of coffea's vector, candidate and NanoAOD classes. No live clash. The only overlaps (`rho`, `m`, `e`, `p`, `deltaR`) are plain output records without a vector behaviour, and `pvec` in `BTA_producer.py` is set on `JetPFCands` records, which are not vectors.
- **Unrelated to the migration**: `Rereco17_94X` and `Winter22Run3` are missing from `campaign_map()` (§4), so `load_SF` fails for both. They are kept on purpose; adding them to `campaign_map()` would make them loadable.

## 10. Top risks (post-introspection — reduced)
1. **numpy 2.x** behavioral changes across ~40 files (`np.float_`/`np.bool8`/`np.NaN` removals, copy-on-write, `np.unique` API) — now the #1 risk since all coffea classes survived. Repo-wide numpy-2 audit. **✅ Resolved:** caused zero failures across all sweeps (§8b), and the alias audit is clean (§9 P3).
2. **JEC internals**: classes survive but `.build()` dropped `lazy_cache`. ~~the pickled factory blob must be regenerated~~ — **superseded by §4**: the pickle path is dormant, so no regeneration is needed. The hand-rolled JES/JER block (`correction.py ~1300-1600`) has since been exercised by every sweep run.
3. **`hist` 2.10 drift** affecting `plot_utils.py` reimplemented `plotratio` / `._storage_type()` internals. `._storage_type()` is replaced by the public `storage_type()`; `plotratio` is still unchecked.
4. **`from_root` default mode='virtual'**: confirm `Runner`+`IterativeExecutor` feeds materializing arrays (virtual/eager), not dask, so `len(events)`/`.to_numpy()` keep working. (Use non-dask executors for first validation.) **✅ Resolved:** eager idioms worked from the first run with both `IterativeExecutor` and `FuturesExecutor` (§8b).

_(Resolved by introspection: column_accumulator/set_accumulator, Rochester txt_converters/rochester_lookup, LumiMask, extractor, BTagScaleFactor, Weights, ProcessorABC — all still present, no longer risks.)_
