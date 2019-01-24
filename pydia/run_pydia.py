#! /usr/bin/env python
from __future__ import print_function

#@author Camilo Jimenez (IAC), Steve Rodney (USC)
#This script is quick+dirty implementation of the Python Difference Imaging
# Analysis (pyDIA) pipeline.
# You can use this script to run pydia from command line, e.g.,
#python run_pyDIA.py --n_parallel=20 --gain=1.9 --readnoise=5.0 --name_pattern=*.fits --pdeg=2 --sdeg=2 --bdeg=2 --datekey=None  --loc_data=liverpool_lens/BELLSJ0747+5055/images_stack --loc_output=liverpool_lens/BELLSJ0747+5055/output_2018-12-02_v1256  --ref_image_list=ref_file.lst --use_GPU=True --min_ref_images=1 --psf_fit_radius=5.0 --iterations=3
#for a complete list of parameters and more, see the pyDIA documentation:
# http://www2.phys.canterbury.ac.nz/%7Emda45/pyDIA/pyDIA-documentation.pdf

#TODO : replace getopt with argparse
#TODO : update documentation

import sys, getopt
import os
import datetime
import time
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.nddata.utils import Cutout2D
import glob

def spin_and_trim(imlist, refimname, trimfrac=0.4, verbose=False):
    """ Rotate images to match the WCS of the reference image (spin) and then
    cut off a fraction of the outer region of the image (trim).
    Returns a list of the trimmed images.
    """
    trimfrac = trimfrac
    if verbose:
        print('Spinning Input Images: ' + str(imlist))
        print("to match WCS of Ref image: " + refimname)
        print("and trimming by : %i pct"%(trimfrac*100))

    trimmed_image_list = []
    hdr_imref = fits.getheader(refimname)
    wcs_imref = WCS(hdr_imref)
    for imname in imlist:
        imname_trimmed = imname.replace('.fits', '_trim.fits')
        if os.path.isfile(imname_trimmed):
            print("%s exists.  Skipping trimming."%imname_trimmed)
        if verbose:
            print("Reprojecting %s to x,y frame of %s " % (
                imname, refimname))
        im = fits.open(imname)
        if imname != refimname:
            array, footprint = reproject_interp(im[0], hdr_imref)
        else:
            # Ref image does not need reprojection
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

        # TODO : handle this with parameters
        # move input image into sub-folder
        imname_moved = os.path.join('UNTRIMMED', os.path.basename(imname))
        if not os.path.exists('UNTRIMMED'):
            os.makedirs('UNTRIMMED')
        os.rename(imname, imname_moved)
        trimmed_image_list.append(imname_trimmed)
    return(trimmed_image_list)


def runpydia(parameters):
    #from pydia import calibration_functions as cal

    params =parameters

    #
    # Import the high-level pipeline routines #
    use_GPU = params.use_GPU
    if use_GPU:
        from pydia import DIA_GPU as DIA
    else:
        from pydia import DIA_CPU as DIA

    if params.ref_image:
        # User has provided a single image to use as the reference image, so
        # write this into the refimagelist and update the params.
        params.ref_image_list = os.path.join(params.loc_data,
                                             'ref_image_list.txt')
        fout = open(params.ref_image_list, 'w')
        print(params.ref_image, file=fout)
        fout.close()

    init=datetime.datetime.now()
    time_int =  time.time()
    log_file= os.path.join(parameters.loc_output, "pydia_run_log.log")
    if not os.path.exists(parameters.loc_output):
        os.makedirs(parameters.loc_output)
    f=open(log_file,"a+")
    f.write(str(init)+" processing "+params.loc_data +" in "+params.loc_output +"\n")
    f.close()

    # first, spin and trim the input images
    #try:
    if True:
        image_list = glob.glob(os.path.join(params.loc_data, params.name_pattern))
        #TODO : accommodate cases where more than one ref image provided
        trimmed_image_list = spin_and_trim(image_list, params.ref_image,
                                           params.trimfrac)
        # TODO : write something useful to the log file
    # except:
    else:
        out=datetime.datetime.now()
        time_out =  time.time()
        f=open(log_file,"a+")
        f.write(str(out)+" error "+str(sys.exc_info()[0])+params.loc_data +" in "+params.loc_output +" "+ str((time_out-time_int)/60)+ " min \n")
        f.close()

    # Next, run PyDIA to make diff images (no photometry!)
    #try:
    if True:
        # after spinning and trimming, we only want to process _trim.fits files
        # TODO: make this better!
        if params.name_pattern.endswith("*.fits"):
            params.name_pattern.replace("*.fits","*_trim.fits")

        if not os.path.exists(params.loc_output):
            DIA.imsub_all_fits(params)

        # TODO : do we need to do the calibrate step?
        # If so, need to include more dependencies, like sklearn
        #cal.calibrate(params.loc_output)
        out=datetime.datetime.now()
        time_out =  time.time()
        f=open(log_file,"a+")
        f.write(str(out)+" Success!  Processed "+params.loc_data +" into "+params.loc_output +" in "+ str((time_out-time_int)/60)+ " min \n")
        f.close()
        # TODO: clean up after successful completion
    #except:
    else:
        out=datetime.datetime.now()
        time_out =  time.time()
        f=open(log_file,"a+")
        f.write(str(out)+" error "+str(sys.exc_info()[0])+params.loc_data +" in "+params.loc_output +" "+ str((time_out-time_int)/60)+ " min \n")
        f.close()


def viewHelp():
    print("For a complete documentation go to http://www2.phys.canterbury.ac.nz/%7Emda45/pyDIA/pyDIA-documentation.pdf")
    print("#Parameter (default value) description")
    print("-h --help help")
    print("-i, --loc_data (images) Images location folder")
    print("-o, --loc_output (images) Output folder, there should no be. This will be create")
    print("--bdeg (0) Degree of spatial variation of the differential bakground")
    print("--ccd_group_size (100) For CCD code version, only reevaluate the PSF for position changes greater than this parameter.")
    print("--cluster_mask_radius (50) Radius of circular region to mask if mask_cluster = True.")
    print("--datekey ('MJD-OBS') Date field in FITS headers")
    print("--detect_threshold (4.0) Sigma threshold for variable object detection.")
    print("--diff_std_threshold (10.0) Threshold for good images to use for variable object detection")
    print("--do_photometry (True) Measure stellar fluxes from the difference images")
    print("--fft_kernel_threshold (3.0)")
    print("--fwhm_mult (6.5) Multiplier to determine kernel size")
    print("--fwhm_section (None) Array of 4 numbers describing the bottom-left and top-right corners of a rectangular section of each image to use for FWHM estimation")
    print("--gain (1.0) Inverse-gain of the CCD (e-/ADU)")
    print("--image_list_file (images) Write image names and dates.")
    print("--iterations (1) Number of kernel iterations")
    print("--kernel_maximum_radius (20.0) Maximum radius for the convolution kernel")
    print("--kernel_minimum_radius (5.0) Minimum radius for the convolution kernel")
    print("--loc_data (.) Absolute or relative path to the directory containing the input images.")
    print("--loc_output (.) Absolute or relative path to the directory to store the output files. Will be created if it doesnt exist.")
    print("--make_difference_images (True) ")
    print("--mask_cluster (False) If True, the attempt to detect and mask the highest concentration of stars in each image.")
    print("--min_ref_images (3) Minimum number of images to combine for the reference.")
    print("--n_parallel (1) Number of parallel CPU parallel processes to run. Typically set this to equal the number of CPU cores available.")
    print("--name_pattern (*.fits) Pattern describing data file names")
    print("--nstamps (200) How many stamps to use")
    print("--pdeg (0) Degree of spatial variation of the kernel photmetric scale")
    print("--pixel_max (50000) Maximum valid pixel value")
    print("--pixel_min (0.0) Minumum valid pixel value")
    print("--pixel_rejection_threshold (3.0) Threshold for masking outlying pixels after each kernel iteration.")
    print("--preconvolve_images (False)")
    print("--preconvolve_FWHM (1.5)")
    print("--psf_fit_radius (3.0) Radius (in pixels) for PSF photometry.")
    print("--psf_profile_type (gaussian) The only option at present.")
    print("--readnoise (1.0) CCD readout noise (e)")
    print("--ref_image (None) Provide a single image to use as the reference image (supersedes ref_include_file).")
    print("--ref_image_list (ref.images) List of images used to create the reference.")
    print("--ref_include_file (None) If this file exists, it should contain a list of images to use for the photometric reference")
    print("--ref_exclude_file (None) If this file exists, it should contain a list of images to exclude from the photometric reference")
    print("--reference_min_seeing (1.0) Minimum FWHM for images to be included in the photometric reference")
    print("--reference_max_roundness (1.3)")
    print("--reference_seeing_factor (1.01) Include images in the photometric reference that have FWHM less than this factor times the lowest FWHM.")
    print("--reference_sky_factor (1.3) Include images in the photometric reference that have backgrounds less than this factor times that of the image with the lowest FWHM.")
    print("--registration_image (None)Use this FITS file as the astrometric reference. Otherwise use the best-seeing image.")
    print("--sdeg (0) Degree of spatial variation of the kernel shape variation")
    print("--sky_degree (0) Degree of spatial variation allowed for the sky background model")
    print("--sky_subtract_mode (percent) If this parameter is set to default, fit a 2D polynomial model for the sky. If this parameter is set to percent, then subtract a constant percentage value from each image.")
    print("--sky_subtract_percent (0.01) Sky percentage to subtract if sky_subtract_mode = percent")
    print("--stamp_edge_distance (40) Minimum distance in pixels from the centre of a stamp to the edge of the detector.")
    print("--stamp_half_width (20) Size of each stamp")
    print("--star_detect_sigma (12) Minimum signal-to-noise for star detection.")
    print("--star_file (None) If provided, use this catalogue of star positions for photometry. The first 3 columns must be star_number, x_position, y_position. Rows starting with # are ignored.")
    print("--star_file_has_magnitudes (False) If true, assume column 4 of star_file contains magnitudes.")
    print("--star_file_is_one_based (True) If true, assume the coordinates in star_file are based on (1,1) being the lower left pixel of the detector.")
    print("--star_file_number_match (10000)")
    print("--star_file_transform_degree 2")
    print("--star_reference_image None")
    print("--subtract_sky (False) Pre-subtract the sky background from each image")
    print("--trimfrac (0.4) Fraction of the image to trim off after spinning to match ref image WCS and before pyDIA processing")
    print("--use_fft_kernel_pixels (False)")
    print("--use_GPU (True) Flag to indicate GPU or CPU use.")
    print("--use_stamps (False) Use small image sections to compute the kernel (as opposed to using the whole image)")
    print("--trim (True) Trim the images to ")

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def isBoolean(value):
    if value == "True" or value == "False":
        try:
            bool(value)
            return True
        except ValueError:
            return False
    else:
        return False

def isInt(value):
    if str(value).find(".") == -1:
        try:
            int(value)
            return True
        except ValueError:
            return False
    else:
        return False
    


def getTypeValue(value):
    val = value
    if isInt(value):
        return int(value)
    elif isBoolean(value):
        return True if value == "True" else False
    elif isfloat(value):
        return float(value)
    elif val == "None":
        return None
    else:
        return value

def main(argv):
    from pydia import data_structures as DS

    try:
        opts, args = getopt.getopt(
            argv, "hvi:o:",
            ["bdeg=", "ccd_group_size=", "cluster_mask_radius=", "datekey=", "detect_threshold=",
            "diff_std_threshold=", "do_photometry=", "fft_kernel_threshold=", "fwhm_mult=",
            "fwhm_section=", "gain=", "image_list_file=", "iterations=",
            "kernel_maximum_radius=", "kernel_minimum_radius=", "loc_data=", "loc_output=",
            "make_difference_images=", "mask_cluster=", "min_ref_images=", "n_parallel=",
            "name_pattern=", "nstamps=", "pdeg=", "pixel_max=", "pixel_min=",
            "pixel_rejection_threshold=", "preconvolve_images=", "preconvolve_FWHM=",
            "psf_fit_radius=", "psf_profile_type=", "readnoise=", "ref_image=",
             "ref_image_list=",
            "ref_include_file=", "ref_exclude_file=", "reference_min_seeing=",
            "reference_max_roundness=", "reference_seeing_factor=", "reference_sky_factor=",
            "registration_image=", "sdeg=", "sky_degree=", "sky_subtract_mode=",
            "sky_subtract_percent=", "stamp_edge_distance=", "stamp_half_width=",
            "star_detect_sigma=", "star_file=", "star_file_has_magnitudes=",
            "star_file_is_one_based=", "star_file_number_match=", "star_file_transform_degree=",
            "star_reference_image=", "subtract_sky=", "trimfrac=",
             "use_fft_kernel_pixels=", "use_GPU=",
            "use_stamps=", "verbose"])
    except getopt.GetoptError:
        print('run_pyDIA.py -i <Imagefolder> -o <Outputfolder>')
        sys.exit(2)

    parameters = DS.Parameters()

    # Set some defaults appropriate for the LCOGT data processing:
    parameters.name_pattern = '*.fits'
    parameters.trimfrac = 0.4
    parameters.min_ref_images = 1
    parameters.use_GPU = False
    parameters.do_photometry = False
    parameters.psf_fit_radius = 5
    parameters.use_stamps = True
    parameters.nstamps = 50
    parameters.subtract_sky = True
    parameters.pixel_min = 0.
    parameters.pixel_max = 50000.

    parameters.loc_input = 'DIA_IN'
    parameters.loc_data = 'DIA_IN'
    parameters.loc_output = 'DIA_OUT'

    for opt, arg in opts:
        if opt in ('-h', "--help"):
            viewHelp()
            sys.exit()
        elif opt in ("-v","--verbose"):
            parameters.verbose = True
        elif opt in ("-i", "--loc_data"):
            parameters.loc_data = arg
        elif opt in ("-o", "--loc_output"):
            parameters.loc_output = arg
        else:
            item = str(opt).replace("--", "")
            val=getTypeValue(arg)
            parameters.__dict__[item] = val

    # If the REF subdir exists, assume it holds the desired ref image
    refdir = os.path.join(parameters.loc_data, 'REF')
    if os.path.isdir(refdir):
        refimlist = glob.glob(os.path.join(refdir,"*.fits"))
        parameters.ref_image = refimlist[0]

    runpydia(parameters)

if __name__ == "__main__":
   main(sys.argv[1:])