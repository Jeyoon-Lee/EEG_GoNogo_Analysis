## Environment

- **MATLAB** with **EEGLAB v2025.1.0**

EEGLAB plugins used:

| Plugin | Version |
| --- | --- |
| ADJUST | v1.1.1 |
| EEG-BIDS | v10.3 |
| FASTER | v1.2.4 |
| ICLabel | v1.7 |
| MFFMatlabIO | v5.0 |
| clean_rawdata | v2.11 |
| dipfit | v5.6 |
| firfilt | v2.8 |

## Scripts

### Pipeline / preprocessing

- **`MADE_revised.m`** — Top-level driver: sets `cfg` (paths, filters, epoching) and
  runs the MADE preprocessing pipeline (filter → ICA → epoching) over the raw dataset.

### ERP analysis

- **`Compute_ERP.m`** — Computes grand-average ERP waveforms and mean amplitudes over a
  frontocentral ROI (FA/error vs Hit/correct), and runs the condition, directionality,
  and group×condition (mixed ANOVA) statistics.
- **`ERP_window_sensitivity.m`** — Loads subject waveforms once and re-measures amplitude
  across several measurement windows to check that condition and group×condition effects
  are consistent regardless of window choice.
- **`ERP_minfa_sensitivity.m`** — Same one-time load, but sweeps the epoch-count screening
  criteria (min_fa / min_go / min_nogo) to see how sample size and statistics change.
- **`plot_erp_ridge.m`** — Ridge plot of subject-averaged ERPs with channels ordered
  anterior→posterior, plotting FA vs Hit to see how a peak (e.g. parietal Pe) varies front-to-back.

### `utils/`

- **`get_MADE_filtered_data.m`** — MADE STEP 1–7: import, channel locations, filtering,
  FASTER bad-channel detection, and E65 reference removal; caches the filtered `.set`.
- **`get_MADE_ica_data.m`** — MADE STEP 8–11: runs ICA and rank-aware ADJUST artifact
  rejection (accepts IC < channels), computing `icaact` manually.
- **`get_MADE_processed_data.m`** — MADE STEP 12–16: response/stimulus-locked epoching with
  condition labels (GoNogo/Accuracy/RT) and fiducial-aware 64-channel interpolation.
- **`Add_block_event.m`** — Tags each `EEG.event` with block / trial-in-block / global-trial
  numbers using `bgin` (trial) and TRSP→bgin gaps (block) as boundaries.
- **`generate_perform_metrics.m`** — Reads block/trial-tagged filtered data and computes
  Go/No-Go performance metrics (Hit / OE / FA / CR, reaction times) into a CSV report.
- **`plot_raw.m`** — Visual inspection of an EEGLAB dataset via FieldTrip `ft_databrowser`,
  with optional channel highlighting.

## Attribution

The EEG preprocessing is based on the **MADE pipeline** (Child Development Lab,
University of Maryland). The `Adjusted_adjust_scripts/` and `appendix_scripts/`
folders and the MADE pipeline functions belong to their original authors and are
redistributed under **GNU GPL v2**. See the original source below:

- **MADE pipeline** — Debnath et al. (2020), *Psychophysiology*, e13580.
  Original repository: https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline

The **ERP analysis scripts** (ERP computation, sensitivity analyses,
plotting) are my own work.

## Modifications

Changes made to the original MADE / ADJUST code:

**`Adjusted_adjust_scripts/`** — made rank-aware so ADJUST works when the number
of ICs is smaller than the number of channels. Only array dimensions were changed; the detection
logic and thresholds are unchanged.
- `beall_horizontal.m`, `beall_blink_detection.m`, `Spatial_Info_eyes.m`:
  `zeros(n,n)` → `zeros(n,nchannels)`, channel loops `1:n` → `1:nchannels`.
- `MARA_extract_time_freq_features.m`: removed per-component plot/JPG saving.

**`get_MADE_filtered_data.m`** (STEP 1–7): minor path/config handling; FASTER
bad-channel detection and E65 reference removal retained.

**`get_MADE_ica_data.m`** (STEP 8–11): rank-aware ADJUST integration — accept
rectangular icaweights (IC < channels), compute `icaact` manually, wrap ADJUST
in try-catch; downsample ICA copy to 250 Hz for speed.

**`get_MADE_processed_data.m`** (STEP 12–16): response/stimulus-locked epoching
with condition labels (GoNogo/Accuracy/RT); removed redundant frontal epoch
rejection (handled by ADJUST); fiducial-aware 64-channel interpolation
(type=='EEG', exclude E65); condition-wise epoch counts.