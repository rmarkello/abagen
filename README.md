# abagen

A small tool to make it easier to work with the Allen Brain Institute [human micorarray expression](http://human.brain-map.org/microarray/search) data.

## Usage

At a minimum, you need an atlas image that you want to use to parcellate the microarray expression data.
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

Once you have an atlas and, if desired, a supplementary CSV file, you can download the AHBA data by calling:

```python
>>> files = abagen.fetch_microarray(donors='all')
```

If you do not specify `donors='all'` microarray expression data from only one donor will be downloaded (instead of all six).
If you have already downloaded the microarray expression from the Allen Brain Institute website, you can set the `data_dir` argument to use those files:

```python
>>> files = abagen.fetch_microarray(data_dir='/path/to/my/download', donors='all')
```

The returned object is a dictionary pointing to the different data files supplied by the Allen Institute.
We can then use those files, along with the `atlas` and `atlas_info`, to generate a region x gene microarray expression array:

```python
>>> expression = abagen.get_expression_data(files, atlas, atlas_info)
```
