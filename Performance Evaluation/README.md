# Wideband Spectrum Sensing: Baseline and Deep Learning Pipelines

This repository contains the implementation of wideband spectrum sensing
using both classical signal processing techniques and deep
learning-based semantic segmentation. The framework enables reproducible
experiments for baseline detection, deep neural network training,
dataset generation, theoretical robustness analysis, and performance
evaluation.

The implementation is developed in MATLAB and relies on standard
communication and deep learning toolboxes.

------------------------------------------------------------------------

## Repository Structure

    ├── BaselineDetection/
    ├── BaselineDetectionMethodPipeline_version2.0/
    ├── Bounds/
    ├── DataCreation/
    ├── DeepLearningDetection/
    ├── PerformanceEvaluation/

Each module is described below.

------------------------------------------------------------------------

## 1. Baseline Detection

This module provides a reference implementation of a classical spectrum
sensing pipeline. It demonstrates spectrogram generation for multiple
wireless technologies and performs stripe-level detection using energy
and cyclostationary features.

### Technologies Supported

-   WLAN (via `wlanNonHTConfig`)
-   LTE (via `lteRMCDL`)
-   NR (via 5G NR Toolbox)
-   RADAR (via Phased Array System Toolbox)

### Important Notes

-   The template scripts use a default sampling rate of **300 MHz**.
-   For full experimental reproduction, the sampling rate should be set
    to **61.44 MHz**, with corresponding center frequencies and
    bandwidth parameters updated accordingly.

### Entry Script

``` matlab
generate_spectrumSensingBaseline.m
```

------------------------------------------------------------------------

## 2. Baseline Detection Method Pipeline (Version 2.0)

This folder contains the structured experimental pipeline used for
large-scale simulations.

### Original Experimental Configuration

-   Number of spectrogram frames: 2000
-   SNR range: \[0, 10, 20, 30, 40\] dB
-   Sampling rate: 61.44 MHz
-   Frame duration: 40 ms

The uploaded template version is configured for demonstration
purposes: - Sampling rate: 300 MHz - Fixed SNR - 180 frames

To reproduce the full experimental setup, these parameters should be
restored.

### Entry Script

``` matlab
genDataset180_and_baseline.m
```

------------------------------------------------------------------------

## 3. Bounds (Theoretical Analysis)

This module derives analytical bounds on the Intersection-over-Union
(IoU) metric under nuisance perturbations such as jitter or minor
geometric distortions. The analysis follows a risk inflation framework
and produces bounds that are generic and applicable beyond spectrum
sensing to broader semantic segmentation tasks.

------------------------------------------------------------------------

## 4. Data Creation

This module generates labeled spectrogram datasets for deep learning
training and evaluation. It supports multiple configurations including:

-   Different SNR levels
-   Variable bandwidths
-   Overlap and non-overlap signal cases
-   Multi-technology coexistence scenarios

All data generation scripts follow the naming convention:

``` matlab
genSpecSegDataset_*.m
```

Users should select the appropriate configuration based on the desired
experiment.

------------------------------------------------------------------------

## 5. Deep Learning Detection

This module implements the semantic segmentation network for wideband
signal identification.

It includes: - Training scripts - Testing scripts - Batch execution
utilities - Model configuration options

The implementation follows MATLAB's Deep Learning Toolbox and supports
segmentation architectures such as DeepLabv3+.

Some helper utilities required for segmentation can be obtained from
MATLAB's official documentation:

https://se.mathworks.com/help/comm/ug/spectrum-sensing-with-deep-learning-to-identify-5g-and-lte-signals.html

------------------------------------------------------------------------

## 6. Performance Evaluation

This module evaluates and compares the classical baseline detector and
the deep learning model.

It supports: - Per-class (one-vs-rest) ROC curves - Macro-averaged ROC
curves - Confusion matrices - AUC computation - F1-score analysis

Before running evaluation scripts: 1. Execute the baseline detection
pipeline. 2. Train and test the deep learning model.

------------------------------------------------------------------------

## Reproducibility Workflow

To reproduce the complete experimental pipeline:

1.  Generate datasets using `DataCreation/` scripts.
2.  Run baseline detection using
    `BaselineDetectionMethodPipeline_version2.0/`.
3.  Train and evaluate the deep learning model in
    `DeepLearningDetection/`.
4.  Generate performance metrics and plots via `PerformanceEvaluation/`.

### Required MATLAB Toolboxes

-   Communications Toolbox
-   LTE Toolbox
-   5G Toolbox
-   WLAN Toolbox
-   Phased Array System Toolbox
-   Deep Learning Toolbox

------------------------------------------------------------------------

## Scope

This repository supports research and experimentation in:

-   Wideband spectrum sensing
-   Multi-technology coexistence (LTE, NR, WLAN, RADAR)
-   Classical baseline detection methods
-   Deep learning-based semantic segmentation
-   ROC/AUC-based evaluation
-   Robustness bounds under nuisance perturbations

The repository is structured to ensure modularity, traceability, and
reproducibility of experiments.
