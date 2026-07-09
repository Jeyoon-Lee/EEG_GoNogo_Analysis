## Attribution

The EEG preprocessing is based on the **MADE pipeline** (Child Development Lab,
University of Maryland). The `Adjusted_adjust_scripts/` and `appendix_scripts/`
folders and the MADE pipeline functions belong to their original authors and are
redistributed under **GNU GPL v2**. See the original source below:

- **MADE pipeline** — Debnath et al. (2020), *Psychophysiology*, e13580.
  Original repository: https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline

The **ERP analysis scripts** (ERN/Pe computation, sensitivity analyses,
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