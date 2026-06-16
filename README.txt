================================================================================
  Analysis pipeline - Groot et al., 2026 
  Netherlands Institute for Neuroscience (NIN)
================================================================================

OVERVIEW
--------
This repository contains the MATLAB analysis pipeline used to reproduce all
figures in:

  Groot, Schall & Westerberg (2026). Visual cortical dynamics supporting predictable attentional capture.

The pipeline loads and preprocesses laminar multi-unit activity (MUA) and
current source density (CSD) recordings from a 15-channel laminar probe in
macaque area V4, then reproduces all main figures.


REQUIREMENTS
------------
- MATLAB R2021b or later
- Signal Processing Toolbox (for filtfilt)
- Statistics and Machine Learning Toolbox (for signrank, kruskalwallis, fitrm)
- Data files (NPY + CSV format) stored on a local or network disk
  (see path configuration in main.m)


REPOSITORY STRUCTURE
--------------------
  main.m          — Step 1: Load and preprocess MUA, CSD, and metadata.
                    Run this first. Produces MUA_datamat, CSD_datamat,
                    metadatamat, and params in the workspace.

  analysis.m      — Step 2: Reproduce all figures.
                    Requires the variables produced by main.m.

  functions/      — Helper functions called by the pipeline:

    Data loading
      matread_npy.m                 read NumPy .npy files into MATLAB
      matread_csv.m                 read CSV metadata files into MATLAB

    Preprocessing
      preprocess_neural_data.m      data loading, baseline correction,
                                    z-scoring, smoothing, and metadata assembly
      baseline_correction.m         normalisation-based baseline correction
      baseline_correction_first.m   first-trial-of-block baseline correction

    Aggregation & statistics
      level_of_analysis.m           trial averaging by session/compartment/channel
      confidence_interval.m         mean ± 95% CI
      running_rs.m                  running Wilcoxon rank-sum test across time
      running_rs_null.m             running rank-sum test against a null distribution
      running_kwt.m                 running Kruskal-Wallis test across time
      permutation_test.m            cluster-based permutation test (Maris & Oostenveld 2007, J. Neuro. Methods)
      pagesLTest.m                  Page's L test for monotonic trends

    Sequence analysis
      extract_target_distractor_sequences.m   extract T-D-D-D trial sequences
      filter_sequences_by_length.m            filter sequences by minimum length

    Visualisation
      plot_significance_bar.m       overlay significance bars on time traces
      plot_chan_data.m              adaptation sequence boxplots + statistics
      plot_chan_diff_session.m      session-demeaned early vs late epoch comparison
      plot_boxplot_lines.m          boxplot with overlaid individual session lines
      SMOOTH_2D.m                   2-D spatial smoothing for CSD heatmaps
      rt_plot.m                     reaction time CDF figure (Figure 1C)
      accuracy_plot.m               accuracy over trials since feature change (Figure 1D)


HOW TO RUN
----------
1. Open MATLAB and set the working directory to this folder (Groot_Westerberg).

2. Edit the path block at the top of main.m to point to your data and output directories:
     disk_dir        = 'D:\path\to\your\data';
     figure_save_dir = 'D:\path\to\output\figures';

3. Run main.m:
     >> run main.m
   This will populate MUA_datamat, CSD_datamat, metadatamat, and params.

4. Run analysis.m:
     >> run analysis.m
   Figures are saved automatically to:
     figures/MUA/PNG/   — rasterised MUA figure exports
     figures/MUA/SVG/   — vector MUA figure exports
     figures/CSD/PNG/   — rasterised CSD figure exports
     figures/CSD/SVG/   — vector CSD figure exports


FIGURE INDEX
------------
  Figure 1C       — Reaction time CDFs (predictable vs unpredictable)
  Figure 1D       — Accuracy over trials since feature change
  Figure 1E       — RT kernel density estimates, quartilized
  Figure 2        — Per-channel laminar MUA traces (target vs distractor)
  Figure 3A       — Target-selection profile by RT quartile, per compartment
  Figure 3B       — Priming effect (predictable - unpredictable) per compartment
  Figure 4B       — Histogram of T-D-D-D sequence counts per session
  Figure 4C       — Neural adaptation MUA traces (T -> D1 -> D2 -> D3 -> D4)
  Figure 4D       — Early/late epoch activity boxplots per laminar compartment
  Figure 5A       — Laminar MUA quartile traces (target)
  Figure 5C       — CSD heatmaps + CSD traces (target)
  Figure 6A       — Laminar MUA quartile traces (distractor)
  Figure 6B       — CSD heatmaps + CSD traces (distractor)


LAMINAR COMPARTMENTS
--------------------
  Channels  1– 5 : Upper
  Channels  6–10 : Middle
  Channels 11–15 : Lower


NOTES
-----
- Sessions 10–29 were recorded in microvolts (uV); the pipeline automatically
  converts CSD values to nA/mm3 by dividing by 1000.
- Smoothing uses a causal-symmetric boxcar filter via filtfilt (zero-phase),
  with a 15-sample window for MUA and 1-sample window for CSD by default.
- Color preference (preferred vs non-preferred feature color) is estimated
  session-wise via a population reliability permutation test (1000 draws of
  100 trials) and stored in metadatamat.pref (1 = pref., 0 = no pref., -1 = non pref.).
  NOTE that color preference is not taken into account in the current analyses.
================================================================================
