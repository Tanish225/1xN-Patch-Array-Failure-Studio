# CST Field Export Processing & Array Analysis

Python-based parsing, post-processing, and visualization pipeline for CST Microwave Studio 3D electric field vector exports (`E field.txt`).

## Overview

This repository provides automated routines to parse Cartesian field vector data exported from CST Studio Suite simulations, calculate complex field magnitudes, and generate spatial distribution maps across planar cross-sections.

## Data Format Specifications

The parser handles CST ASCII exports structured in a multi-plane planar grid ($x \in [-65.39, 65.39]\text{ mm}$, $y \in [-65.39, 65.39]\text{ mm}$) spanning discrete Z-planes:

| Column Index | Field Name | Description | Units |
| --- | --- | --- | --- |
| **1–3** | `x`, `y`, `z` | Spatial Cartesian coordinates | mm
| **4–6** | `ExRe`, `EyRe`, `EzRe` | Real parts of E-field components ($\text{Re}\{E_x\}, \text{Re}\{E_y\}, \text{Re}\{E_z\}$) | V/m
| **7–9** | `ExIm`, `EyIm`, `EzIm` | Imaginary parts of E-field components ($\text{Im}\{E_x\}, \text{Im}\{E_y\}, \text{Im}\{E_z\}$) | V/m


## Formulations

The total complex electric field vector $\vec{E}(\vec{r})$ is reconstructed from real and imaginary vector components:

$$\vec{E}(\vec{r}) = (\text{ExRe} + j\cdot\text{ExIm})\hat{a}_x + (\text{EyRe} + j\cdot\text{EyIm})\hat{a}_y + (\text{EzRe} + j\cdot\text{EzIm})\hat{a}_z$$



The peak field magnitude $\vert{}\vec{E}\vert{}$ at each grid point is computed using:

$$\vert{}\vec{E}\vert{} = \sqrt{ (\text{ExRe}^2 + \text{ExIm}^2) + (\text{EyRe}^2 + \text{EyIm}^2) + (\text{EzRe}^2 + \text{EzIm}^2) }$$



## Project Structure

```text
├── arrayfailurestudio.m
├── LICENSE                # Open-source license (MIT/Apache 2.0)
└── README.md              # Documentation

```

## Features

* **Multi-Plane Data Ingestion:** Automated chunking of 3D exports into sliced Z-plane grids.


* **Complex Vector Computation:** Vectorized NumPy calculation of component amplitudes and overall magnitude $\vert{}\vec{E}\vert{}$.


* **Visualization Suite:**
* 2D heatmap color meshes of spatial E-field distributions per Z-plane slice.


* 1D centerline intensity profiles along principal axes ($x=0$ / $y=0$).





## Getting Started

### Installation

```bash
git clone https://github.com/Tanish225/cst-field-array-analysis.git
cd cst-field-array-analysis
pip install -r requirements.txt

```
