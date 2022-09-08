# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017

Function that . . . 
            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


def Background_Corrected(Random_Sample_Name, Field_2, Match_Radius, Initial_Background, Magnitude_Column, N_Bins, Range_Min, Range_Max, N_Sample, N_Original_First_Field, Match_Rescale):
    from Topcat_Match_self import skymatch_self
    from Topcat_Match_remove import skymatch_remove
    from Topcat_Column_Remove import Topcat_Column_Remove
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import math

    Original_Match_Radius = Match_Radius
    New_Match_Radius      = Match_Radius * Match_Rescale
    if Match_Rescale != 1:
        print('\n\nWARNING: Match_Rescale == ', Match_Rescale, '\n\n')
    
    #Performing the self match. Removing the central, original sources
    skymatch_self(
        Random_Sample_Name,
        ["RA_MOCK", "DEC_MOCK"],
        Field_2,
        ["RA", "DEC"],
        New_Match_Radius, 'log/Random_Sample_Match_Temp.fits')

    skymatch_remove(
        'log/Random_Sample_Match_Temp.fits',
        ["RA", "DEC"],
        Random_Sample_Name,
        ["RA_MOCK", "DEC_MOCK"],
        0.01, 'log/Random_Sample_Match.fits')

       

    #Removing the original columns to remove confusion
    Topcat_Column_Remove('log/Random_Sample_Match.fits', "*_1", 'log/Random_Sample_Match.fits')

    Corrected_Background_Sources = fits.open('log/Random_Sample_Match.fits')[1].data
    Corrected_Magnitudes         = Corrected_Background_Sources[Magnitude_Column]
    print('\nCorrected_Background_Length          == ', len(Corrected_Magnitudes))
    Corrected_Magnitudes         = Corrected_Magnitudes[np.logical_not(np.isnan(Corrected_Magnitudes))]
    print('Corrected_Background_Length (No NaN) == ', len(Corrected_Magnitudes), '\n')    
    Histogram_Corrected_F, Histogram_Corrected_Edges = np.histogram(Corrected_Magnitudes, bins=N_Bins, range=(Range_Min,Range_Max))

    Initial_Background_Edges = Initial_Background[:,0]
    Initial_Background_F     = Initial_Background[:,1]

    print('\nDoing the background rescale:')
    print('N_Original_First_Field    == ', N_Original_First_Field)
    print('N_Sample                  == ', N_Sample)
    print('Match_Rescale             == ', Match_Rescale)
    print('Background Rescale Factor == ', 1  / (Match_Rescale**2) *N_Original_First_Field / N_Sample  )
    print('\n')
 
    Corrected_Background_Edges = Histogram_Corrected_Edges[:-1]  #Removing the last unnecessary edge
    Corrected_Background_F     = Histogram_Corrected_F / (Match_Rescale**2) *N_Original_First_Field / N_Sample
    Corrected_Background       = np.column_stack((  np.array(Corrected_Background_Edges[:]), np.array(Corrected_Background_F[:])  ))

   

    #Creating the function output
    col1 = fits.Column(name='Left-Bin_Edge',        format='D', array=Initial_Background_Edges)
    col2 = fits.Column(name='Initial_Background',   format='D', array=Initial_Background_F)
    col3 = fits.Column(name='Corrected_Background', format='D', array=Corrected_Background_F)

    cols = fits.ColDefs([col1, col2, col3])
    
    Fits_File = fits.BinTableHDU.from_columns(cols)
    Fits_File.writeto('log/Histogram_Data_Background.fits', overwrite=True)




    return {'Corrected_Background': Corrected_Background}


    
