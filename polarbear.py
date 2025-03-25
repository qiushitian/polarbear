#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simplistic reduction scripts for ARCTIC on the APO 3.5m telescope.

Author: Qiushi (Chris) Tian
Last updated: 2025-03-25
"""
import numpy as np
import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty


def trim_overscan(
        data,
        n=1024,
        x_overscan_start=2,
        y_overscan_start=0,
        x_overscan_center=50,
        y_overscan_center=2
):
    """
    Output a data array with overscan pixels truncated.

    :param data: The image *data* to be trimmed (not a HDU list or HDU).
    :param n: Quadrant pixel size.
        Do not change if running on ARCTIC data. Default: 1024.
    :param x_overscan_start: First non-overscan pixel in the *x*-direction.
        Do not change if running on ARCTIC data. Default: 2.
    :param y_overscan_start: First non-overscan pixel in the *y*-direction.
        Do not change if running on ARCTIC data. Default: 0.
    :param x_overscan_center: Number of overscan pixels in the center of the frame in the *x*-direction.
        Do not change if running on ARCTIC data. Default: 50.
    :param y_overscan_center: Number of overscan pixels in the center of the frame in the *y*-direction.
        Do not change if running on ARCTIC data. Default: 2.
    :return: Trimmed data.
    """
    # old unused code
    # np.append(data[: n, :], data[n + y_overscan_center :, :], axis=0)

    ll = data[
         y_overscan_start: y_overscan_start + n,
         x_overscan_start: x_overscan_start + n
         ]
    lr = data[
         y_overscan_start: y_overscan_start + n,
         x_overscan_start + n + x_overscan_center: x_overscan_start + n + x_overscan_center + n
         ]
    ul = data[
         y_overscan_start + n + y_overscan_center: y_overscan_start + n + y_overscan_center + n,
         x_overscan_start: x_overscan_start + n
         ]
    ur = data[
         y_overscan_start + n + y_overscan_center: y_overscan_start + n + y_overscan_center + n,
         x_overscan_start + n + x_overscan_center: x_overscan_start + n + x_overscan_center + n
         ]

    trimmed = np.zeros((n * 2, n * 2), dtype=float)
    trimmed[:n, :n] = ll
    trimmed[:n, n:] = lr
    trimmed[n:, :n] = ul
    trimmed[n:, n:] = ur

    return trimmed


def make_master_bias(bias_list, output_path):
    """
    Write a master bias file from raw bias files.

    :param bias_list: A list of paths to raw bias files.
        This can be something like: `[f'bias/Bias.00{_}.fits' for _ in range(24, 34 + 1)]`.
    :param output_path: Write path for the output master bias file.
    """
    ccdproc.combine(
        bias_list, output_path,
        clip_extrema=True, nlow=1, nhigh=1,
        sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
        sigma_clip_func='median', signma_clip_dev_func='mad_std', unit='adu'
    )


def flatfield(sci_path: str, bias_path, flat_path: str, write_path: str, overwrite=False):
    """
    Performs flat field correction.

    :param sci_path: Path to the untrimmed science file.
    :param bias_path: Path to the untrimmed master bias file.
    :param flat_path: Path to the untrimmed flat file.
    :param write_path: Write path for the output file.
    :param overwrite: Whether to overwrite existing output file. Default: False.
    """
    # bias
    bias_hdul = fits.open(bias_path)
    bias_ccdd = CCDData(
        bias_hdul['PRIMARY'].data,
        unit='adu',
        uncertainty=StdDevUncertainty(bias_hdul['UNCERT'].data, unit='adu')
    )

    # sci
    sci_ccdd = CCDData(fits.getdata(sci_path), unit='adu')
    sci_ccdd = ccdproc.subtract_bias(sci_ccdd, bias_ccdd)

    # flat
    flat_ccdd = CCDData(fits.getdata(flat_path), unit='adu')
    flat_ccdd = ccdproc.subtract_bias(flat_ccdd, bias_ccdd)

    corrected_ccdd = ccdproc.flat_correct(sci_ccdd, flat_ccdd)
    trimmed = trim_overscan(corrected_ccdd.data)
    fits.writeto(write_path, trimmed, overwrite=overwrite)


if __name__ == '__main__':
    # REPLACE THESE EXAMPLES!
    make_master_bias(
        [f'path/to/biases/Bias.00{_}.fits' for _ in range(24, 34 + 1)], 'write/path/master_bias.fits'
    )

    flatfield(
        'path/to/sci_filter_1.fits', 'path/to/master_bias.fits', 'path/to/flat_filter.fits', 'write/path/reduced_1.fits'
    )
    flatfield(
        'path/to/sci_filter_2.fits', 'path/to/master_bias.fits', 'path/to/flat_filter.fits', 'write/path/reduced_2.fits'
    )
    flatfield(
        'path/to/sci_filter_3.fits', 'path/to/master_bias.fits', 'path/to/flat_filter.fits', 'write/path/reduced_3.fits'
    )
