# -*- coding: utf-8 -*-

import os
import itertools

import nibabel as nib
import numpy as np
from scipy.ndimage import center_of_mass
from scipy.spatial.distance import cdist

from .transforms import ijk_to_xyz


def leftify_atlas(atlas):
    """
    Zeroes out all ROIs in the right hemisphere of `atlas`

    Assumes that positive X values indicate the right hemisphere (e.g., RAS+
    orientation) and that the X-origin is in the middle of the brain

    Parameters
    ----------
    atlas : str or niimg-like
        Filepath to or in-memory loaded image

    Returns
    -------
    atlas : niimg-like
        Loaded image with right hemisphere zeroed out
    """

    atlas = check_img(atlas)

    # get ijk corresponding to zero-point
    i, j, k = nib.affines.apply_affine(np.linalg.inv(atlas.affine), [0, 0, 0])

    # zero out all positive voxels; img is RAS+ so positive = right hemisphere
    data = np.array(atlas.dataobj, copy=True)
    data[int(i):] = 0

    return atlas.__class__(data, atlas.affine, header=atlas.header)


def check_img(img):
    """
    Very basic checker that loads `img`` and ensures it's 3D/int

    Parameters
    --------
    img : str or niimg-like
        Filepath to or in-memory loaded image

    Returns
    -------
    img : niimg-like
        Loaded 3D/int image
    """

    if isinstance(img, (str, os.PathLike)) and os.path.exists(img):
        img = nib.load(img)
    elif not isinstance(img, nib.spatialimages.SpatialImage):
        raise TypeError('Provided image must be an existing filepath or a '
                        'pre-loaded niimg-like object')

    # ensure 3D or squeezable to 3D
    if len(img.shape) == 4 and img.shape[3] == 1:
        data = np.asarray(img.dataobj)
        affine = img.affine
        img = img.__class__(data[:, :, :, 0], affine, header=img.header)
    elif len(img.shape) != 3:
        raise ValueError('Provided image must be 3D')

    # check if atlas data is int or castable to int
    # if image is arrayproxy convert it to an array for speed-up
    data = np.asarray(img.dataobj)
    cast = nib.is_proxy(img.dataobj)
    if img.header.get_data_dtype().kind not in ['i', 'u']:
        idata = data.astype('int32')
        cast = np.allclose(idata, data)
        data = idata
        if not cast:
            raise ValueError('Provided image should have integer values or '
                             'be safely castable to int without data loss')
    if cast:
        img = img.__class__(data, img.affine, header=img.header)
        img.header.set_data_dtype(np.int32)

    return img


def get_unique_labels(label_image):
    """
    Returns all possible ROI labels from ``label_image``

    Parameters
    ----------
    label_image : niimg-like object
        ROI image, where each ROI is identified with a unique integer ID

    Returns
    -------
    labels : np.ndarray
        Integer labels of all ROIS found within ``label_image``
    """

    label_image = check_img(label_image)
    return np.trim_zeros(np.unique(label_image.dataobj)).astype(int)


def get_centroids(image, labels=None, image_space=False):
    """
    Finds centroids of ``labels`` in ``label_image``

    Parameters
    ----------
    label_image : niimg-like object
        3D image containing integer label at each point
    labels : array_like, optional
        List of values containing labels of which to find centroids. Default:
        all possible labels
    image_space : bool, optional
        Whether to return xyz (image space) coordinates for centroids based
        on transformation in ``label_image.affine``. Default: False

    Returns
    -------
    centroids : (N, 3) np.ndarray
        Coordinates of centroids for ROIs in input data
    """

    image = check_img(image)
    data = np.asarray(image.dataobj)

    # if no labels of interest provided, get all possible labels
    if labels is None:
        labels = np.trim_zeros(np.unique(data))

    # get centroids for all possible labels
    centroids = np.row_stack(center_of_mass(data, labels=data, index=labels))

    # return xyz if desired; otherwise, ijk
    if image_space:
        centroids = ijk_to_xyz(centroids, image.affine)

    return centroids


def closest_centroid(coord, centroids, return_dist=False):
    """
    Returns index of `centroids` closest to `coord`

    Computes Euclidean distances between `coord` and each of `centroids`,
    returning index of closest centroid

    Parameters
    ----------
    coord : (1, 3) array_like
        Coordinates of sample
    centroids : (N, 3) array_like
        Centroids of parcels (in same space as `coord`)
    return_dist : bool, optional
        Whether to also return distance of closest centroid

    Returns
    -------
    closest : int
        Index of closest centroid in `centroids` to `coord`
    distance : float
        Distance of closest centroid in `centroids` to `coord`. Only returned
        if `return_dist=True`
    """

    distances = np.squeeze(cdist(np.atleast_2d(coord), centroids))
    closest = distances.argmin(axis=0)

    if return_dist:
        return closest, distances[closest]

    return closest


def expand_roi(coords, dilation=1, return_array=True):
    """
    Expands coordinates ``coords`` to include neighboring coordinates

    Computes all possible coordinates of distance ``dilation`` from ``coords``.
    Returns a generator (``itertools.product``) by default, but can return an
    array if desired (if ``return_array=True``).

    Parameters
    ----------
    coords : (1, 3) array_like
        List of ijk values for coordinate in 3D space
    dilation : int, optional
        How many neighboring voxels to expand around `coords`. Default: 1
    return_array : bool, optional
        Whether to return array instead of generator. Default: True

    Returns
    -------
    coords : (27, 3) generator
        Coordinates of expanded ROI
    """

    def expand(x, d=1):
        return np.arange((x - d), (x + d + 1), dtype=int)

    # return all combinations of coordinates
    coords = np.squeeze(coords)
    gen = itertools.product(*[expand(n, d=dilation) for n in coords])

    # coerce to array if desired
    if return_array:
        return np.asarray(list(gen))

    return gen
