[![GitHub code licence is MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)]

# Matlab source code to calculate Coefficient of Variation given HiTC Hi-C matrix

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)

## Overview

Hi-C and 5C data are noisy, especially long range interactions. To filter these data, we implemented a new method based on coefficient of variation to flag noisy interactions

## Installation

```bash
git clone git@github.com:zhanyinx/Coefficient_Variation.git
```

## Usage

Modify the parameters in the main script coefficient_variation.m according to your dataset.
Create an excel table:
Sample	WT/Mut	Type of Mut	Size	'Name'	chrom	chromStart	chromEnd

| Sample   | WT/Mut | Type of Mut | Size | 'Name' | chrom	| chromStart | chromEnd  |
| -------  | ------ | ----------- | ---- | ------ | ----- | ---------- | --------  |
| B20_18C6 | Mut    | Inversion   | 70kb | 70kb.  | chrX  | 100625999  | 100701714 |
| .......  | ...... | ........... | .... | ...... | ..... | .......... | ........  |
