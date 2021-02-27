# -*- coding: utf-8 -*-
"""
Functions for transforming between coordinate systems
"""

from pathlib import Path

import nibabel as nib
import numpy as np

from .datasets import fetch_freesurfer, fetch_raw_mri

MNI152TO305 = np.array([[1.0022, 0.0071, -0.0177, 0.0528],
                        [-0.0146, 0.9990, 0.0027, -1.5519],
                        [0.0129, 0.0094, 1.0027, -1.2012]])
MNI305TO152 = np.array([[0.9975, -0.0073, 0.0176, -0.0429],
                        [0.0146, 1.0009, -0.0024, 1.5496],
                        [-0.0130, -0.0093, 0.9971, 1.1840]])


def ijk_to_fsnative(ijk, donor, data_dir=None):
    """
    Converts provided `ijk` coordinates to RAS fsnative space of `donor`

    Fsnative space is returned with :func:`abagen.fetch_freesurfer`

    Parameters
    ----------
    ijk : (N, 3) array_like
        IJK coordinates to be transformed to RAS surface space
    donor : str
        Which donor `coords` belongs to
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. Default: $HOME/
        abagen-data

    Returns
    -------
    ras : (N, 3) np.ndarray
        RAS coordinates
    """

    ijk = np.asarray(ijk)

    # load orig.mgz volume from freesurfer and get torig
    orig = fetch_freesurfer(donors=donor, data_dir=data_dir, verbose=0)[donor]
    mgz = nib.load(Path(orig) / 'mri' / 'orig.mgz')
    torig = mgz.header.get_vox2ras_tkr()

    # load raw mri from AHBA and get affine
    raw_mri = fetch_raw_mri(donors=donor, data_dir=data_dir, verbose=0)[donor]
    raw_affine = nib.load(raw_mri['t1w']).affine

    # convert ijk to xyz (in raw mri space) then to ijk (in orig.mgz space)
    ijk = np.c_[
        xyz_to_ijk(ijk_to_xyz(raw_affine, ijk), mgz.affine),
        np.ones(len(ijk))
    ]

    # convert ijk (in orig.mgz) to RAS (in fsnative)
    return ijk @ torig[:-1].T


def mni152_to_fsaverage(xyz):
    """
    Converts MNI152 `xyz` coordinates to fsaverage RAS (MNI305) coordinates

    Parameters
    ----------
    xyz : (N, 3) array_like
        MNI152 coordinates

    Returns
    -------
    ras : (N, 3) np.ndarray
        fsaverage (MNI305) coordinates
    """

    xyz = np.asarray(xyz)
    return np.c_[xyz, np.ones(len(xyz))] @ MNI152TO305.T


def fsaverage_to_mni152(xyz):
    """
    Converts fsaverage RAS (MNI305) `xyz` coordinates to MNI152 coordinates

    Parameters
    ----------
    xyz : (N, 3) array_like
        fsaverage (MNI305) coordinates

    Returns
    -------
    ras : (N, 3) np.ndarray
        MNI152 coordinates
    """

    xyz = np.asarray(xyz)
    return np.c_[xyz, np.ones(len(xyz))] @ MNI305TO152.T


def ijk_to_xyz(coords, affine):
    """
    Converts `coords` in cartesian space to `affine` space

    Parameters
    ----------
    coords : (N, 3) array_like
        Cartesian (ijk) coordinate values
    affine : (4, 4) array_like
        Affine matrix containing displacement + boundary

    Returns
    -------
    xyz : (N, 3) np.ndarray
        Provided ``coords`` in ``affine`` space
    """

    return nib.affines.apply_affine(affine, coords)


def xyz_to_ijk(coords, affine):
    """
    Converts `coords` in `affine` space to cartesian space

    Parameters
    ----------
    coords : (N, 3) array_like
        Image coordinate (xyz) values
    affine : (4, 4) array_like
        Affine matrix containing displacement + boundary

    Returns
    ------
    ijk : (N, 3) numpy.ndarray
        Provided `coords` in cartesian space
    """

    return nib.affines.apply_affine(np.linalg.inv(affine), coords).astype(int)
