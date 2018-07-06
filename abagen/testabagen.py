import abagen

in_dir = ('/home/gshafiei/data/Projects/TEMPLATES_rois_desikan/')

scale = [33, 60, 125, 250, 500]
nnodes = [83, 129, 234, 463, 1015]

s = scale[ii]
my_img = (in_dir + 'ICBM152_nlin_sym_09a/ROIv_scale%s_RAS_ICBM152.nii.gz' % s)

roi_img = my_img  # each ROI should have unique integer ID

# download microarray data for all available donors
files = abagen.fetch_microarray(donors='all')

# get expression data for each ROI in ``roi_img``, averaging across donors
expression = abagen.get_expression_data(files, roi_img, metric='mean')
