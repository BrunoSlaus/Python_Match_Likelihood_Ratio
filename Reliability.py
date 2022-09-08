# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017

Function that, for a given matching-radius, calculates the reliability
(see Will Sutherland and Will Saunders, 1992, MNRAS, 259, 413-420) for
each matched sources. 

The function will create the MATCHED.fits file that contains the 
matched sources and their reliability. It will return the total
number of matched sources and the inputed radius.
            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


def Reliability(Q):
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    from scipy.interpolate import interp1d
    
    Matched_LR = fits.open('Output/Matched_LR.fits')[1].data
    LR         = Matched_LR['Likelihood_Ratio']
    GroupID    = Matched_LR['GroupID'] 
    GroupSize  = Matched_LR['GroupSize']

    REL_Array = np.array([])
    for i in np.arange(0, len(LR)): 
        if GroupSize[i]>1:
            Group = Matched_LR[Matched_LR['GroupID']==GroupID[i]]
            Sum = np.sum(Group['Likelihood_Ratio'])
        else:
            Sum = LR[i]

        REL = LR[i] / (Sum + (1-Q))
        REL_Array = np.append(REL_Array, REL)
        


    REL_Column  = fits.Column(name='Reliability', array=REL_Array, format='D')
    REL_Columns = fits.ColDefs([REL_Column])
    Data_Columns = Matched_LR.columns
    Merged_Columns = Data_Columns + REL_Columns



    Fits_File = fits.BinTableHDU.from_columns(Merged_Columns)
    Fits_File.writeto('Output/Matched_REL.fits', overwrite=True)
    print('Output saved into Output/Matched_REL.fits file.')
 
    return 








