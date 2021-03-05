# -*- coding: utf-8 -*-
"""
Functions for transforming between coordinate systems
"""

from pathlib import Path

import nibabel as nib
import numpy as np

from .datasets import fetch_freesurfer

MNI152TO305 = np.array([[1.0022, 0.0071, -0.0177, 0.0528],
                        [-0.0146, 0.9990, 0.0027, -1.5519],
                        [0.0129, 0.0094, 1.0027, -1.2012],
                        [0.0000, 0.0000, 0.0000, 1.000]])
MNI305TO152 = np.linalg.inv(MNI152TO305)


def _get_fs_affine_torig(donor, data_dir=None):
    """
    Gets FreeSurfer orig.mgz affine and torig matrices for `donor`

    Parameters
    ----------
    donor : str
        Which donor to get data from
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. This function
        will fetch the FreeSurfer directory of `donor` to get the necessary
        information. Default: $HOME/abagen-data

    Returns
    -------
    affine : (4, 4) np.ndarray
        Affine matrix of orig.mgz
    torig : (4, 4) np.ndarray
        Torig matrix of orig.mgz
    """

    # load orig.mgz volume from freesurfer and get torig
    orig = fetch_freesurfer(donors=donor, data_dir=data_dir, verbose=0)[donor]
    mgz = nib.load(Path(orig) / 'mri' / 'orig.mgz')
    torig = mgz.header.get_vox2ras_tkr()

    return mgz.affine, torig


def xyz_to_fsnative(xyz, donor, data_dir=None):
    """
    Converts provided `xyz` image coordinates to RAS fsnative space of `donor`

    Parameters
    ----------
    xyz : (N, 3) array_like
        XYZ coordinates (in `donor` image space) to be transformed to RAS
        surface space
    donor : str
        Which donor `xyz` coordinates belong to
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. This function
        will fetch the FreeSurfer directory of `donor` to get the necessary
        information for converting the coordinates. Default: $HOME/abagen-data

    Returns
    -------
    ras : (N, 3) np.ndarray
        RAS coordinates
    """

    xyz = np.atleast_2d(xyz)
    affine, torig = _get_fs_affine_torig(donor, data_dir=data_dir)
    return ijk_to_xyz(xyz_to_ijk(xyz, affine, floor=False), torig)


def fsnative_to_xyz(fsnative, donor, data_dir=None):
    """
    Converts provided `fsnative` RAS coordinates to xyz image space for `donor`

    Parameters
    ----------
    fsnative : (N, 3) array_like
        RAS coordinates (in `donor` fsnative space) to be transformed to xyz
        image space
    donor : str
        Which donor `fsnative` coordinates belong to
    data_dir : str, optional
        Directory where data should be downloaded and unpacked. This function
        will fetch the FreeSurfer directory of `donor` to get the necessary
        information for converting the coordinates. Default: $HOME/abagen-data

    Returns
    -------
    xyz : (N, 3) np.ndarray
        XYZ coordinates
    """

    fsnative = np.atleast_2d(fsnative)
    affine, torig = _get_fs_affine_torig(donor, data_dir=data_dir)
    return ijk_to_xyz(xyz_to_ijk(fsnative, torig, floor=False), affine)


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

    xyz = np.atleast_2d(xyz)
    return ijk_to_xyz(xyz, MNI152TO305)


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

    xyz = np.atleast_2d(xyz)
    return ijk_to_xyz(xyz, MNI305TO152)


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

    coords, affine = np.atleast_2d(coords), np.asarray(affine)
    return nib.affines.apply_affine(affine, coords)


def xyz_to_ijk(coords, affine, floor=True):
    """
    Converts `coords` in `affine` space to cartesian space

    Parameters
    ----------
    coords : (N, 3) array_like
        Image coordinate (xyz) values
    affine : (4, 4) array_like
        Affine matrix containing displacement + boundary
    floor : bool, optional
        Whether to round down converted `ijk` coordinates to int. Default: True

    Returns
    ------
    ijk : (N, 3) numpy.ndarray
        Provided `coords` in cartesian space
    """

    coords, affine = np.atleast_2d(coords), np.asarray(affine)
    ijk = nib.affines.apply_affine(np.linalg.inv(affine), coords)
    if floor:
        ijk = np.asarray(np.floor(ijk), dtype=int)
    return ijk
