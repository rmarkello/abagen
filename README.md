# abagen

A small tool to make it easier to work with the Allen Brain Institute [human micorarray expression](http://human.brain-map.org/microarray/search) data.

## Usage

```python
import abagen

roi_img = 'my_roi_image.nii.gz'  # each ROI should have unique integer ID

# download microarray data for all available donors
files = abagen.fetch_microarray(subjects='all')

# get expression data for each ROI in ``roi_img``, averaging across donors
expression = abagen.get_expression_data(files, roi_img, metric='mean')
```
