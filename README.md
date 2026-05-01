# Data_analysis

MATLAB scripts for preprocessing sorted electrophysiology data, extracting trial-wise stimulus metadata, preparing DLAG-ready inputs, training DLAG models, and analyzing inferred latent structure.

This repository is organized around a practical end-to-end workflow:

1. align spikes/events to trials
2. extract stimulus information from Expo XML
3. compute run-level and condition-level unit metrics
4. bin and normalize neural activity
5. build `model_data_allruns.mat`
6. train DLAG models
7. analyze inferred latents and condition effects

---

## Overview

This repo is a personal working codebase for multi-probe neural data analysis centered on **DLAG (Delayed Latents Across Groups)**.

It combines:
- session-level preprocessing from SpikeGLX / CatGT / Kilosort outputs
- trial-wise stimulus parsing from Expo XML
- unit-level quality and response summaries
- construction of DLAG-ready MATLAB data structures
- DLAG training in pooled or condition-specific modes
- post hoc analysis of latent reproducibility, latent categories, shared variance, and condition dependence

The code is currently organized as **script-based workflows**, not as a packaged toolbox. In most cases, you should open a script, edit the user parameters near the top, and run it for your session.

---

## Repository structure

```text
data_seg/
  sorting_check/
    get_spikt_evt_mrkidx_seg.m
    unit_statistic_by_condition.m
  generate_lfp_seg_trial.m
  processing_to_count_and_fr.m

readstim_info/
  stiminfo_pertrial.m

dlag_mainheng/
  prepar_data/
    model_data_prepar.m
  train/
    train_dlag.m
    train_dlag_by_condition.m
  analysis/
    Latents_compare.m
    Anova_latents_for_all_conds_used_dlag.m
    plot_dlag_results.m
