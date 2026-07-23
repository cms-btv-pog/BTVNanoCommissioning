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
| **Execution** | `runner.py` (only) | `run_uproot_job` + `iterative/futures/parsl/dask_executor` | **Rewrite** to `Runner` + executor classes (§3) |
| **JEC/JER build** | `utils/correction.py` (`JME_shifts`, ~1300–1669), `helpers/func.py` (factory) | `CorrectedJetsFactory.build(..., lazy_cache=events.caches[0])` | Drop `lazy_cache`; regen factory pickle (§4) |
| **Accumulators** | `helpers/func.py` (`dump_lumi`) | `column_accumulator`, `set_accumulator` | ✅ No change (both survive) |
| **Histogram fills** | `utils/histogramming/histogrammer.py` | eager flatten + numpy | Minimal; fix latent `if <ak> > <ak>` bug (§6) |
| **Env/pins** | `setup.cfg`, `test_env.yml` | coffea 0.7.31 | DONE |
| **Outliers** | `BTA_producer.py`, `BTA_ttbar_producer.py`, `sf_ttdilep_kin.py` | extra `np.array(ak…)`, BDT `to_numpy` | Validate under eager |

Everything else (corrections `extractor`/`LumiMask`/`Weights`/`BTagScaleFactor`, plotting scripts, condor submitters) works with minor/no changes — confirmed by introspection (§7).

## 2. Environment (done)
- `setup.cfg`: `coffea[dask,dask-awkward,xrootd]==2026.7.0`, `correctionlib>=2.6.0`, `python_requires=>=3.10`.
- `test_env.yml`: python 3.10–3.12, coffea 2026.7.0 via pip, `parsl>=2024.12.09`, dropped `setuptools<=70.1.1` cap.
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
✅ **No change required** — `column_accumulator`, `set_accumulator`, and `accumulate` all survive in 2026.7.0 (§7). Optional cleanup: remove the dead `@property accumulator`/`self._accumulator` (never assigned) from all 18 processors when touching them.

## 6. Histogramming — `utils/histogramming/histogrammer.py`
- Fills already eager: `h.fill(temp_syst, flatten(genflavor), flatten(pruned_ev.SelJet[histname]), weight=flatten(...))`, `syst=np.full(len(weight), syst)`. Works under eager/virtual NanoEvents.
- `hist.Hist.fill` (hist 2.10) accepts the flattened awkward + numpy weights (API stable).
- **Latent bug to fix regardless**: `histogrammer.py:400-416` uses `pt if <ak_array> > <ak_array> else ...` — a Python `if` on a non-scalar awkward array. Rarely reached today; fix with `ak.where`.

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
- Env naming: canonical name stays `btv_coffea` (CI `activate-environment: btv_coffea`; hardcoded in runner.py:613, condor_lxplus/submitter.py:231/235, scripts/fetch.py). Keep the migrated production env named `btv_coffea`; `btv_coffea_2026` is a dev sandbox only.
- Postprocessing scripts (`scripts/plotdataMC.py`, `make_template.py`, `comparison.py`, `plotSysts.py`, …) load `.coffea` via `coffea.util.load`; output object shape unchanged in eager mode → expected to keep working (watch `hist` version drift in `plot_utils.py`).
- CI workflows (`.github/workflows/*_workflow.yml`) build `test_env.yml` and run one workflow each — these become the end-to-end migration gate.

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

## 9. Phasing / PR breakdown
- **P0 (done)**: pins + env + `import coffea` smoke test.
- **P1**: `runner.py` execution rewrite; get `example` workflow end-to-end on 1 file with `IterativeExecutor` (eager). Green = events load + histograms fill + `.coffea` saved.
- **P2**: JEC/JER — drop `events.caches[0]`, regenerate factory pickle, fix `.build()` call; get `JME_shifts` (nominal + up/down) working.
- **P3**: sweep the other 17 processors; fix eager outliers (`BTA_producer`, `BTA_ttbar_producer`, `sf_ttdilep_kin` BDT `to_numpy`); numpy-2 audit (`np.float_`, copy semantics).
- **P4**: scale-out executors (dask/parsl/condor) + CI green across all `*_workflow.yml`.
- **P5**: plotting/template scripts + docs.

## 10. Top risks (post-introspection — reduced)
1. **numpy 2.x** behavioral changes across ~40 files (`np.float_`/`np.bool8`/`np.NaN` removals, copy-on-write, `np.unique` API) — now the #1 risk since all coffea classes survived. Repo-wide numpy-2 audit.
2. **JEC internals**: classes survive but `.build()` dropped `lazy_cache`, and the pickled factory blob must be regenerated under coffea 2026 / awkward 2.11. The hand-rolled JES/JER block (`correction.py ~1300-1600`) needs runtime validation.
3. **`hist` 2.10 drift** affecting `plot_utils.py` reimplemented `plotratio` / `._storage_type()` internals.
4. **`from_root` default mode='virtual'**: confirm `Runner`+`IterativeExecutor` feeds materializing arrays (virtual/eager), not dask, so `len(events)`/`.to_numpy()` keep working. (Use non-dask executors for first validation.)

_(Resolved by introspection: column_accumulator/set_accumulator, Rochester txt_converters/rochester_lookup, LumiMask, extractor, BTagScaleFactor, Weights, ProcessorABC — all still present, no longer risks.)_
