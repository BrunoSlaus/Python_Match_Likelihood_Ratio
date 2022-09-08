# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017

Function that, 
            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


def Random_Sample(Input_Field_Name, N_Sample, Magnitude_Distribution, Magnitude_Column, Match_Radius, Move_Rescale):
    from Topcat_Subset import Topcat_Subset
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    from astropy.table import vstack, Table
    import astropy.units as u
    from subprocess import call



    if Move_Rescale != 1:
        print('\n\nWARNING: Move_Rescale == ', Move_Rescale, '\n\n')
    
    Edges      = Magnitude_Distribution[:,0]
    Edges_Diff = np.diff(Edges)[0]
    Edges_Last = Edges[-1] + Edges_Diff
    Edges      = np.append(Edges, Edges_Last)
    
    Normalisation = np.sum(Magnitude_Distribution[:,1])

    print('\nStarting the creation of the Random subsample.')
    print('If the Full-Counterpart field is large, this might take a while.\n')

    
    #Binning the Full-Counterpart data
    Sample_Array = []
    for i in range(0,len(Edges)-1):  
        print('Progress: Bin '+ str(i+1) + ' / ' + str(len(Edges)-1) )
        Min_Edge = Edges[i]
        Max_Edge = Edges[i+1]

        Topcat_Subset(Input_Field_Name, Magnitude_Column, Min_Edge, Max_Edge, 'log/Full_Bin_' + str(i+1) +'.fits')
        Full_Bin = fits.open('log/Full_Bin_' + str(i+1) +'.fits')[1].data
        call('rm -f log/Full_Bin_' + str(i+1) +'.fits', shell=True)  #REMOVING THE BIN ONCE IT WAS USED

        N_Sample_Bin = Magnitude_Distribution[:,1][i] / Normalisation * N_Sample
        N_Sample_Bin = int(N_Sample_Bin)
        if len(Full_Bin)>0 and N_Sample_Bin>0:    #Checking if it is empty
            Full_Bin = np.random.permutation(Full_Bin)
            Sample_Bin = Full_Bin[0:N_Sample_Bin]
            #print('Sample_Array', Sample_Array)
            #print('N_Sample_Bin', N_Sample_Bin)
            #print('Sample_Bin', Sample_Bin['NUMBER'])
            Sample_Array.append(Sample_Bin)
    Sample_Array_Output = np.concatenate(Sample_Array)
    #print(Sample_Array_Output['NUMBER'])

    Fits_File_Unmoved = fits.BinTableHDU.from_columns(Sample_Array_Output)
    Fits_File_Unmoved.writeto('log/Sample_Array_Unmoved.fits', overwrite=True)



    #Performing the random move of each mock source
    N_Length_Sample = len(Sample_Array_Output['RA']) #Very similar to N_Sample
    Theta = np.random.random((N_Length_Sample, 1)) * math.pi * 2
    r     = np.sqrt(np.random.random((N_Length_Sample, 1)) )* Match_Radius * Move_Rescale
    Ra_Move, Dec_Move = [r*np.cos(Theta), r*np.sin(Theta)]
    Ra_Move  = Ra_Move  * u.arcsec
    Dec_Move = Dec_Move * u.arcsec
    Ra_Move  = Ra_Move.to(u.deg).value
    Dec_Move = Dec_Move.to(u.deg).value
    #print("Sample_Array_Output['RA']",Sample_Array_Output['RA'])
    #print('Length', len(Sample_Array_Output['RA']))
    #print(Ra_Move[:,0])
    #print(len(Ra_Move))
    Sample_Array_Output['RA']  = Sample_Array_Output['RA']  + Ra_Move[:,0]
    Sample_Array_Output['DEC'] = Sample_Array_Output['DEC'] + Dec_Move[:,0]

    col1 = fits.Column(name='NUMBER_MOCK', format='D', array=Sample_Array_Output['NUMBER'])
    col2 = fits.Column(name='RA_MOCK',  format='D'   , array=Sample_Array_Output['RA'])
    col3 = fits.Column(name='DEC_MOCK', format='D'   , array=Sample_Array_Output['DEC'])
    col4 = fits.Column(name='MAG_MOCK', format='D'   , array=Sample_Array_Output[Magnitude_Column])

    cols = fits.ColDefs([col1, col2, col3, col4])
    
    Fits_File = fits.BinTableHDU.from_columns(cols)
    Fits_File.writeto('log/Sample_Array.fits', overwrite=True)



    return {'Sample_Array':Sample_Array_Output, 'Random_Sample_Name':'log/Sample_Array.fits' }    






