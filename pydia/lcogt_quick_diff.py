#! /usr/bin/env python
from __future__ import print_function

import sys
import os
import argparse
import datetime
import time
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.nddata.utils import Cutout2D
import glob
import data_structures as DS
import io_functions as IO
import image_functions as IM
import fnmatch
import numpy as np

# TODO : let user define CPU vs GPU
import DIA_CPU as DIA
# import DIA_GPU as DIA


def set_default_parameters():

    parameters = DS.Parameters()

    # Set some defaults appropriate for the LCOGT data processing:
    parameters.name_pattern = '*.fits'
    parameters.trimfrac = 0.4
    parameters.min_ref_images = 1
    parameters.reference_seeing_factor = 1.1
    parameters.use_GPU = False
    parameters.do_photometry = False
    parameters.psf_fit_radius = 5
    parameters.use_stamps = False
    parameters.nstamps = 50
    parameters.subtract_sky = True
    parameters.pixel_min = 0.
    parameters.pixel_max = 50000.

    parameters.reference_seeing_factor = 1.3
    parameters.loc_input = 'DIA_IN'
    parameters.loc_data = 'DIA_IN'
    parameters.loc_trim = 'DIA_TRIM'
    parameters.loc_output = 'DIA_OUT'
    parameters.wcs_ref_image = None
    parameters.ref_image = 'ref.fits'
    parameters.ref_image_list = 'ref_image_list.txt'
    parameters.verbose = False

    parameters.dospintrim = False
    parameters.dorefim = False
    parameters.dodiffs = False

    return(parameters)


def copy_wcs(wcs_source_file, wcs_target_file):
    """ Copy the WCS header keywords from the source file into the
    target file.
    """
    hdr_src = fits.getheader(wcs_source_file)
    wcs_src = WCS(hdr_src)

    im = fits.open(wcs_target_file, 'update')
    im[0].header.update(wcs_src.to_header())
    im.flush(output_verify='fix+warn')
    im.close()
    return


def spin_and_trim(imlist, wcsrefimname, trimfrac=0.4, trimdir='DIA_TRIM',
                  verbose=False):
    """ Rotate images to match the WCS of the WCS reference image (spin) and
    then cut off a fraction of the outer region of the image (trim).
    Returns a list of the trimmed images.
    """
    if verbose:
        print('Spinning Input Images: ' + str(imlist))
        print("to match WCS Ref image: " + wcsrefimname)
        print("and trimming by : %i pct"%(trimfrac*100))

    if not os.path.exists(trimdir):
        os.makedirs(trimdir)
    trimmed_image_list = []
    hdr_imref = fits.getheader(wcsrefimname)
    wcs_imref = WCS(hdr_imref)
    for imname in list(imlist) + [wcsrefimname]:
        imname_trimmed = os.path.join(
            trimdir, os.path.basename(imname).replace('.fits', '_trim.fits'))
        if os.path.isfile(imname_trimmed):
            print("%s exists.  Skipping trimming."%imname_trimmed, file=sys.stderr)
            continue
        if verbose:
            print("Reprojecting %s to x,y frame of %s " % (
                imname, wcsrefimname), file=sys.stderr)
        im = fits.open(imname)
        if imname != wcsrefimname:
            array, footprint = reproject_interp(im[0], hdr_imref)
        else:
            # WCS Ref image does not need reprojection
            array = im[0].data

        # Trim off the outer trimfrac to reduce chances of NaN errors due to
        # empty array segments after rotation (also speeds up DIA processing)
        arraysize = array.shape
        cutout = Cutout2D(array, position=[arraysize[1]/2., arraysize[0]/2.],
                          size=[round(arraysize[1]*(1-trimfrac)),
                                round(arraysize[0]*(1-trimfrac))],
                          wcs=wcs_imref)

        # save reprojected-and-trimmed image:
        im[0].data = cutout.data
        im[0].header.update(cutout.wcs.to_header())
        im.writeto(imname_trimmed, output_verify='fix+warn')
        im.close()

        trimmed_image_list.append(imname_trimmed)
    return(trimmed_image_list)


def get_observation_list(filenamelist, params):
    """ From a given list of filenames, make a list of images as
    pyDIA Observation objects.
    If user has specified a full path to the input image files,
    use that full path.  If not, then check if each input image
    file exists as a trimmed  version, and use that.  If not, fall
    back is to use the original untrimmed image.
    """
    observationslist = []
    for filename in filenamelist:
        if filename.endswith('_trim.fits'):
            ftrimname = os.path.basename(filename)
        else:
            ftrimname = os.path.basename(filename).replace('.fits', '_trim.fits')
        trimfile = os.path.join(params.loc_trim, ftrimname)
        imfile = os.path.join(params.loc_data, os.path.basename(filename))
        for f in [filename, trimfile, imfile]:
            if os.path.exists(f):
                g = DS.Observation(f, params)
                observationslist.append(g)
                break
    return(observationslist)


def register_images(files, params):
    """Register the images in the given files list.
    Registered images are placed in the loc_output dir with prefix r_
    """
    # Have we specified a registration template?
    if params.registration_image:
        reg = DS.Observation(params.registration_image, params)
    else:
        reg = DS.EmptyBase()
        reg.fw = 999.0
        for f in files:
            if (f.fw < reg.fw) and (f.fw > 1.2):
                reg = f
    print('Registration image:', reg.name)

    # Register images
    for f in files:
        if f == reg:
            f.image = f.data
            rf = params.loc_output + os.path.sep + 'r_' + f.name
            IO.write_image(f.image, rf)
        else:
            f.register(reg, params)
            # delete image arrays to save memory
            del f.image
            del f.mask
            del f.inv_variance
        del reg.data
        del reg.image
        del reg.mask
        del reg.inv_variance
    return(files)


def make_ref_image(params):
    """Make a reference image if it doesn't exist.
    Requires a list of input images to combine, in the file given by
    params.ref_image_list (default "ref_image_list.txt").  Assumes these
    images are already registered and trimmed, and can be found in the
    params.loc_trim directory.

    Returns the ref image as a  DIA data_structures Observation object
    """
    refimname = params.ref_image
    refimpath = os.path.join(params.loc_output, refimname)
    if os.path.exists(refimpath):
        refim = DS.Observation(refimpath, params)
        print("Ref image %s exists. Not clobbering."%refimpath)
        return(refim)

    # read in a list of input images for creating the ref image
    if not os.path.exists(params.ref_image_list):
        print("Missing ref image list file %s"%params.ref_image_list)
        return(None)
    fin = open(params.ref_image_list, 'r')
    ref_input_filenames = [f.strip() for f in fin.readlines()]
    fin.close()

    ref_image_list = get_observation_list(ref_input_filenames, params)
    registered_ref_image_list = register_images(ref_image_list, params)

    # Make the photometric reference image if we don't have it.
    # Find stamp positions if required.
    stamp_positions = DIA.make_reference(registered_ref_image_list, params,
                                         reference_image=refimname)
    refim = DS.Observation(refimpath, params)
    # Register the newly made reference image #
    #TODO: ref image registration is superfluous? maybe a symlink would work?
    registered_refimlist = register_images([refim], params)

    # TODO : we need to add a header and WCS?
    copy_wcs(ref_image_list[0].fullname, refimpath)
    copy_wcs(ref_image_list[0].fullname, registered_refimlist[0].fullname)
    return(refim)


def make_diff_images(filenamelist, refim, params):
    """ make a diff image for each file in filenamelist: file - refim.

    filenamelist :  list of filenames for 'target' images
    refim : reference image, either a filename or a DIA Observation object
    params : DIA parameters object
    """
    star_group_boundaries = None
    detector_mean_positions_x = None
    detector_mean_positions_y = None
    star_unsort_index = None
    star_positions = None
    stamp_positions = None
    sky = 0.0

    if isinstance(refim, str) and os.path.exists(refim):
        refim = DS.Observation(refim, params)

    # TODO: investigate what is really being done here:
    #  Apply saturation mask and boxcar blurring to reference image
    mask, _ = IO.read_fits_file(
            params.loc_output + os.path.sep + 'mask_' + refim.name)
    refim.mask = mask
    pm = params.pixel_max
    params.pixel_max *= 0.9
    refim.mask *= IM.compute_saturated_pixel_mask(refim.image, 4, params)
    params.pixel_max = pm
    refim.blur = IM.boxcar_blur(refim.image)
    if params.mask_cluster:
        refim.mask *= IM.mask_cluster(refim.image, refim.mask, params)

    # For each given filename, get a pyDIA observation object
    image_list = get_observation_list(filenamelist, params)

    # Register the images, using the ref image as the registration template,
    # unless the user has specified otherwise
    if not params.registration_image:
        params.registration_image = refim.fullname
    registered_image_list = register_images(image_list, params)

    # make diff images: im - ref
    for im in registered_image_list:
        result = DIA.difference_image(
            refim, im, params,
            stamp_positions=stamp_positions,
            psf_image=params.loc_output + os.path.sep + 'psf.fits',
            star_positions=star_positions,
            star_group_boundaries=star_group_boundaries,
            detector_mean_positions_x=detector_mean_positions_x,
            detector_mean_positions_y=detector_mean_positions_y)
        del im.image
        del im.mask
        del im.inv_variance

        hdr = fits.getheader(im.fullname)
        #  TODO : use astropy fits to propagate header with WCS from parent image
        # Save output images to files
        if isinstance(result.diff, np.ndarray):
            IO.write_image(result.diff,
                           params.loc_output + os.path.sep + 'd_' + im.name,
                           header=hdr)
            #IO.write_image(result.model,
            #               params.loc_output + os.path.sep + 'm_' + f.name)
            #IO.write_image(result.norm,
            #               params.loc_output + os.path.sep + 'n_' + f.name)
            #IO.write_image(result.mask,
            #               params.loc_output + os.path.sep + 'z_' + f.name)


def do_everything(params=None):
    if params is None:
        params = set_default_parameters()

    # Create the output directory if it doesn't exist
    if not (os.path.exists(params.loc_output)):
        os.mkdir(params.loc_output)

    # The degree of spatial shape changes has to be at least as
    # high as the degree of spatial photometric scale
    if (params.sdeg < params.pdeg):
        print('Increasing params.sdeg to ', params.pdeg)
        params.sdeg = params.pdeg

    # Print out the parameters for this run.
    print('Parameters:')
    for par in dir(params):
        print(par, getattr(params, par))
    print('Determining list of input images')

    input_file_basenames = os.listdir(params.loc_input)
    input_file_paths = []
    for f in input_file_basenames:
        if fnmatch.fnmatch(f, params.name_pattern):
            input_file_paths.append(os.path.join(params.loc_input,f))

    if params.dospintrim:
        print("Spinning and trimming input images")
        if not params.wcs_ref_image:
            refimname = params.ref_image
            refimpath = os.path.join(params.loc_output, refimname)
            if os.path.exists(refimpath):
                params.wcs_ref_image = refimpath
        if not params.wcs_ref_image:
            params.wcs_ref_image = input_file_paths[0]
        spin_and_trim(input_file_paths, params.wcs_ref_image, params.trimfrac)

    if params.dorefim:
        print("Making reference image")
        refim = make_ref_image(params)
    else:
        refimpath = os.path.join(params.loc_output, params.ref_image)
        refim = DS.Observation(refimpath, params)

    trim_file_basenames = os.listdir(params.loc_trim)
    trim_file_paths = []
    # read in a list of input images for creating the ref image
    ref_input_filenames = []
    if os.path.exists(params.ref_image_list):
        fin = open(params.ref_image_list, 'r')
        ref_input_filenames = [os.path.basename(f.strip())
                               for f in fin.readlines()]
        fin.close()

    for f in trim_file_basenames:
        if fnmatch.fnmatch(f, params.name_pattern):
            if f not in ref_input_filenames:
                trim_file_paths.append(os.path.join(params.loc_trim,f))
    print("Trimmed files to subtract: %s"%str(trim_file_paths))
    if params.dodiffs:
        print("Making diff images")
        make_diff_images(trim_file_paths, refim, params)
    return


def main():
    parser = argparse.ArgumentParser(
        description='Run the PyDIA difference image construction pipeline')

    # No required positional arguments
    # parser.add_argument('required_arg', help='Required positional argument.')

    # Optional Arguments
    parser.add_argument('--dospintrim', action='store_true',
                        help="Do the spinning and trimming step.")
    parser.add_argument('--dorefim', action='store_true',
                        help="Do the building a ref image step.")
    parser.add_argument('--dodiffs', action='store_true',
                        help="Do the making diff images step.")
    parser.add_argument('--doall', action='store_true',
                        help="Do all the processing steps.")

    parser.add_argument('--indir', type=str, default='DIA_IN',
                        help='Directory with input images.')
    parser.add_argument('--trimdir', type=str, default='DIA_TRIM',
                        help='Directory for trimmed images.')
    #parser.add_argument('--refdir', type=str, default='DIA_REF',
    #                    help='Directory with
    parser.add_argument('--outdir', type=str, default='DIA_OUT',
                        help='Directory for output diff images, etc.')

    parser.add_argument('--wcsrefim', type=str, default='',
                        help='WCS reference image. Default: first input image.')
    parser.add_argument('--refim', type=str, default='ref.fits',
                        help='Subtraction reference image.')
    parser.add_argument('--refimlist', type=str, default='ref_image_list.txt',
                        help='File with list of images to combine'
                        ' to make the subtraction ref image')

    parser.add_argument('--pattern', type=str, default='*.fits',
                        help='Filename pattern to select files to process.')
    parser.add_argument('-v', dest='verbose', action='count', default=0,
                        help='Turn verbosity up (use -v,-vv,-vvv, etc.)')

    argv = parser.parse_args()
    if argv.doall:
        argv.dospintrim = True
        argv.dorefim = True
        argv.dodiffs = True

    params = set_default_parameters()
    params.dospintrim = argv.dospintrim
    params.dorefim = argv.dorefim
    params.dodiffs = argv.dodiffs
    params.loc_input = argv.indir
    params.loc_trim = argv.trimdir
    params.loc_output = argv.outdir
    params.wcs_ref_image = argv.wcsrefim
    params.ref_image= argv.refim
    params.ref_image_list = argv.refimlist
    params.verbose = argv.verbose

    if not np.any([argv.dospintrim, argv.dorefim, argv.dodiffs]):
        print("You haven't specified any step to do. \n"
              "Use --dospintrim, --dorefim, --dodiffs or --doall.\n"
              "To see all options, use --help or -h.")
    else:
        do_everything(params)


if __name__ == '__main__':
    main()
