# abagen

This package provides a Python interface for working with the [Allen Human Brain Atlas Institute](http://human.brain-map.org/microarray/search) (AHBA) microarray expression data.

[![Build Status](https://travis-ci.org/rmarkello/abagen.svg?branch=master)](https://travis-ci.org/rmarkello/abagen)
[![Codecov](https://codecov.io/gh/rmarkello/abagen/branch/master/graph/badge.svg)](https://codecov.io/gh/rmarkello/abagen)

## Table of Contents

If you know where you're going, feel free to jump ahead:

* [Installation](#installation-and-setup)
* [Purpose](#purpose)   
    * [Overview](#overview)   
    * [Development](#development)   
* [Example usage](#usage)
* [How to get involved](#how-to-get-involved)
* [Credits](#credits)

## Installation and setup

This package requires Python >= 3.5. 
Assuming you have the correct version of Python installed, you can install `abagen` by opening a terminal and running the following:

```bash
git clone https://github.com/rmarkello/abagen.git
cd abagen
python setup.py install
```

There are plans (hopes?) to get this set up on PyPi for an easier installation process, but that is a long-term goal!

## Purpose

### Overview

In 2013, the Allen Brain Institute released a dataset containing microarray expression data from human brain tissue, sampled post-mortem from six donors. 
This dataset has offered an unprecedented opportunity to examine the genetic underpinnings of the human brain, and has already yielded novel insight into e.g., [adolescent brain development](http://www.pnas.org/content/113/32/9105.long) and [functional brain organization](http://science.sciencemag.org/content/348/6240/1241.long).

In order to be effectively used in most analyses, the AHBA microarray expression data needs to be, at a minimum, (1) collapsed into regions of interest (e.g., parcels or networks), and (2) combined across donors.
While this may potentially seem trivial, there are numerous analytic choices in these steps that can dramatically influence the resulting data and any downstream analyses.
Indeed, Arnatkevičiūte et al., 2018 provided a thorough treatment of this in a [recent manuscript](https://www.biorxiv.org/content/early/2018/07/30/380089), demonstrating how the techniques and code used to the prepare the raw AHBA data have varied widely across published reports.

The current Python package, `abagen`, aims to provide a fully reproducible pipeline for processing and preparing the AHBA microarray expression data for analysis.

### Development

This package has been largely developed in the spare time of a single graduate student ([`@rmarkello`](https://github.com/rmarkello))&mdash;with help from some incredible [contributors](https://github.com/rmarkello/abagen/graphs/contributors)).
While it would be :sparkles: amazing :sparkles: if anyone else finds it helpful, given the limited time constraints of graduate school, the current package is not currently accepting requests for new features.
However, if you might be interested in contributing yourself, take a look at [how to get involved](#how-to-get-involved)!

## Usage

At a minimum, you will need an atlas image that you want to use to parcellate the microarray expression data.
The supplied atlas should be a Nifti image where regions or parcels are denoted by unique integer IDs.

```python
>>> import abagen
>>> atlas = 'atlas.nii.gz'
```

You can also supply a CSV file with information about the atlas.
This CSV file will constrain the matching of microarray expression samples to anatomical regions.

If supplied, the CSV _must_ have the following columns:

  1. `id`: ID corresponding to integer label in `atlas` image
  2. `hemisphere`: L/R hemispheric designation
  3. `structure`: broad structural class (i.e., 'cortex', 'subcortex', or 'cerebellum')

For example, a valid CSV might look like this:

```python
>>> atlas_info = 'atlas.csv'
>>> pd.read_csv(atlas_info)
   id                    label hemisphere structure
0   1  lateralorbitofrontal_rh          R    cortex
1   2         parsorbitalis_rh          R    cortex
2   3           frontalpole_rh          R    cortex
3   4   medialorbitofrontal_rh          R    cortex
4   5      parstriangularis_rh          R    cortex
...
```

Notice that extra columns (i.e., `label`) are okay, as long as the three required columns as present.

Once you have an atlas and, if desired, a supplementary CSV file, you can download the AHBA data by calling:

```python
>>> files = abagen.fetch_microarray(donors='all')
```

If you do not specify `donors='all'` microarray expression data from only one donor will be downloaded.
If you have already downloaded the microarray expression from the Allen Brain Institute website, you can set the `data_dir` argument to use those files:

```python
>>> files = abagen.fetch_microarray(data_dir='/path/to/my/download', donors='all')
```

The returned object is a dictionary pointing to the different data files supplied by the Allen Institute.
We can then use those files, along with the `atlas` and `atlas_info`, to generate a region x gene expression array:

```python
>>> expression = abagen.get_expression_data(files, atlas, atlas_info)
```

Unfortunately, due to how samples were collected from the donor brains, it is possible that some regions in the atlas may not be represented by any expression data. If you require a full matrix with expression data for _every_ region, you can specify the following:

```python
>>> expression = abagen.get_expression_data(files, atlas, atlas_info, exact=False)
```

By default, `abagen` `will attempt to be as precise as possible in matching microarray samples with brain regions.
Specifying `exact=False` will, at the cost of this precision, ensure that every brain region is matched to _at least_ one sample.
You can investigate other options for modifying how the `expression` array is generated in the documentation by typing `help(abagen.get_expression_data)`.

## How to get involved

We're thrilled to welcome new contributors! 
If you're interesting in getting involved, you should start by reading our [contributing guidelines](CONTRIBUTING.md) and [code of conduct](CODE_OF_CONDUCT.md).
Once you're done with that, you can take a look at our [issues](https://github.com/rmarkello/abagen/issues) to see if there's anything you might like to work on. 

If you've found a bug, are experiencing a problem, or have a question, create a new issue with some information about it!

## Credit

While this package was initially created in early 2018, many of the current functions in the project were inspired by the workflow laid out in [Arnatkevičiūte et al., 2018](https://www.biorxiv.org/content/early/2018/07/30/380089).
As such, if you use this code it would be good to (1) provide a link back to this repository with the version of the code used, and (2) cite their paper:

> Arnatkeviciute, A., Fulcher, B. D., & Fornito, A. (2018). A practical guide to linking brain-wide gene expression and neuroimaging data. bioRxiv, 380089.
