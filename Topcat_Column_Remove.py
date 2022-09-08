#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
import os
from astropy.io import fits
#########################################################################

def Topcat_Column_Remove(file,
             Columns,
             output=None,
             java="java -jar",
             topcat_loc=r"/Users/Astro/Documents/Software/",
             topcat="topcat-full.jar",
             ext=1, header=False
             ):
    if output is None:
        out = 'tempcat.fits'
    else:
        out = output

    topcat = " ".join([java, topcat_loc+topcat, "-stilts"])
    string = " ".join([topcat,
                       "tpipe",
                       "in=%s" % file,
                       "out=%s" % out,
                       """cmd='delcols "%s"'  """ % Columns ])


    print(string)                   
    os.system(string)
    data = fits.getdata(out, ext=ext)
    if output is None:
        os.remove(out)
    if header:
        head = fits.getheader(string, ext=ext)
        return data, head
    else:
        return data

#########################################################################
# Example usage
"""
skymatch(
    "/Users/kresimirtisanic/topcattest/detections_at_3_and_1.4_GHz.fits",
    "Ra", 31, 32, 
    "/Users/kresimirtisanic/topcattest/test2.fits")
"""
#########################################################################
# End of code
