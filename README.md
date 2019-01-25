# pyDIA
pyDIA is a modular python package for performing star detection, difference imaging and photometry on astronomical images.

If you use this code, please cite

<a href="https://doi.org/10.5281/zenodo.268049"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.268049.svg" alt="DOI"></a>


This fork of the pyDIA package (by S.Rodney) was initially designed to use for building difference images of data from the 'lenslowz' project, an observing program on the LCOGT network, led by PI Sukanya Chakrabarti (RIT).

Notable changes include: 

1. there is a new script called run_pydia.py, which is a combination of the command-line setup that Camilo built, plus the extra image registration and trimming steps that I did.

2. I have updated the way the script runs, so that it is now designed to work "out of the box" from a directory that is set up as follows:
  DIA_IN/  : contains the images you want to process
  DIA_IN/REF :  contains a single image that will be the reference image.  All other images get shifted to match the WCS of this one, and it is used as the reference for every diff image.

3.  If you set up your directory that way, then you don't need to feed in any other parameters to get the script to run using the 'default' parameters that I've defined (suitable for most LCOGT data sets).   But you can adjust any of the pydia parameters by passing in command-line arguments.   Use "run_pydia.py --help" to see all the choices and check the pydia documentation for details.

4.   I *think* that you can also import this into another higher level script, so your data monitoring program can use this to run pydia automatically.    For example:

```
from pydia import run_pydia as rp
params = rp.set_default_parameters()
# modify parameters as desired
rp.make_diff_images(params)
```

Note that you'll need to be a little careful that you 
have set your current working directory cleanly, before 
invoking the ``set_default_parameters()`` and 
``make_diff_images()``  functions, since they have 
adopted the pydia convention of assuming the user 
is calling the script from the directory where the
 data are.