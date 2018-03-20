# -*- coding: utf-8 -*-

import itertools
import numpy as np
import scipy.ndimage as snd
import scipy.spatial as ss


def get_unique_labels(label_image):
    """
    Returns all possible parcel labels from `label_image`

    Parameters
    ----------
    label_image : niimg-like object
        Parcel image, where each parcel should be identified with a unique
        integer ID

    Returns
    -------
    labels : np.ndarray
        Integer labels of all parcels found within ``label_image``
    """

    return np.trim_zeros(np.unique(label_image.get_data())).astype(int)


def get_centroids(label_image, labels_of_interest=None, image_space=False):
    """
    Finds centroids of ``labels_of_interest`` in ``label_image``

    Parameters
    ----------
    label_image : niimg-like object
        3D image containing integer label at each point
    labels_of_interest : array_like, optional
        List of values containing labels of which to find centroids. Default:
        all possible labels
    image_space : bool, optional
        Whether to return x,y,z (image space) coordinates for centroids,
        based on transformation in ``label_image.affine``. Default: False

    Returns
    -------
    centroids : (3 x N) np.ndarray
        Coordinates of centroids for ROIs in input data
    """

    # if no labels of interest provided, get all possible labels
    if labels_of_interest is None:
        labels_of_interest = get_unique_labels(label_image)

    image_data = label_image.get_data()
    centroids = []

    # iterate through all the provided labels
    for label in labels_of_interest:
        # create blank array with only locations = `label` set to 1
        temp = np.zeros_like(image_data)
        temp[image_data == label] = 1
        # find centroid for this `label`
        centroids.append(snd.measurements.center_of_mass(temp))

    centroids = np.row_stack(centroids).T

    # return x,y,z if desired; otherwise, i,j,k
    if image_space:
        return ijk_to_xyz(centroids, label_image.affine[:-1])
    else:
        return centroids


def ijk_xyz_input_check(coords):
    """
    Confirms inputs to `ijk_to_xyz()` and `xyz_to_ijk()` are in proper format

    Parameters
    ----------
    coords : array_like

    Returns
    -------
    coords : (3 x N) np.ndarray
    """

    coords = np.atleast_2d(coords)

    if coords.shape[0] != 3:
        coords = coords.T
    if coords.shape[0] != 3:
        raise ValueError("Input coordinates must be of shape (3,N).")

    return coords


def ijk_to_xyz(coords, affine):
    """
    Converts voxel ``coords`` in cartesian space to ``affine`` space

    Parameters
    ----------
    coords : (3 x N) array_like
        i, j, k values
    affine : (3 x 4) array_like
        Affine matrix containing displacement + boundary

    Returns
    -------
    coords : (3 x N) np.ndarray
        Provided ``coords`` in ``affine`` space
    """

    coords = ijk_xyz_input_check(coords)

    return (affine[:, :-1] @ coords) + affine[:, [-1]]


def xyz_to_ijk(coords, affine):
    """
    Converts voxel `coords` in `affine` space to cartesian space

    Parameters
    ----------
    coords : (3 x N) array_like
        x, y, z values
    affine : (3 x 4) array_like
        Affine matrix containing displacement + boundary

    Returns
    -------
    coords : (3 x N) np.ndarray
        Provided `coords` in cartesian space
    """

    coords = ijk_xyz_input_check(coords)

    return np.linalg.solve(affine[:, :-1], coords - affine[:, [-1]])


def expand_roi(coords, dilation=1, return_array=False):
    """
    Expands coordinates ``coords`` to include neighboring coordinates

    Computes all possible coordinates of distance ``dilation`` from ``coords``.
    Returns a generator (``itertools.product``) by default, but can return an
    array if desired (if ``return_array=True``).

    Parameters
    ----------
    coords : (3 x 1) array_like
        List of i, j, k values for coordinate in 3D space
    dilation : int, optional
        How many neighboring voxels to expand around `coords`. Default: 1
    return_array : bool, optional
        Whether to return generator (default) or array. Default: False

    Returns
    -------
    coords : (27, 3) generator
        Coordinates of expanded ROI
    """

    def to_three(x, d=1):
        return np.arange((x - d), (x + d + 1), dtype='int')

    # return all combinations of coordinates
    gen = itertools.product(*[to_three(n, d=dilation) for n in coords])
    # coerce to array if desired
    if return_array:
        return np.asarray(list(gen))
    else:
        return gen


def closest_centroid(coords, centroids):
    """
    Returns label for parcel with centroid closest to ``coords``

    Computes Euclidean distances between ``coords`` and each centroid in
    ``centroids``, returning index of closest centroid

    Parameters
    ----------
    coords : (1 x 3) array_like
        Coordinates of sample
    centroids : (N x 3) array_like
        Centroids of parcels (in same space as `coords`)

    Returns
    -------
    label : int
        Index of closest centroid in ``centroids``
    """

    distances = np.squeeze(ss.distance.cdist(np.atleast_2d(coords), centroids))

    return distances.argmin()
