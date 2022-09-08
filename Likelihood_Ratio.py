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


def Likelihood_Ratio(Matched_Field_Name, Magnitude_Column, Final_Real, Final_Background, ID_1, Error_Ra_1_Column, Error_Dec_1_Column, Error_Ra_2_Column, Error_Dec_2_Column, N_Original_First_Field, Match_Radius, Q):
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    from scipy.interpolate import interp1d
    
    print('N_Original_First_Field == ', N_Original_First_Field)

    Data = fits.open(Matched_Field_Name)[1].data

    ID         = Data[ID_1]
    Magnitudes = Data[Magnitude_Column]
    Separation = Data['Separation']
    Error_Ra_1  = Data[Error_Ra_1_Column]
    Error_Dec_1 = Data[Error_Dec_1_Column]

    Fixed_Positional_Error = 0.3
    if Error_Ra_2_Column == 'fixed':
        print('\nWARNING: Error_Ra_2 set to fixed value.')
        print('Fixed_Positional_Error = ', Fixed_Positional_Error)
        Error_Ra_2  = Fixed_Positional_Error
    else:    
        Error_Ra_2  = Data[Error_Ra_2_Column]

    if Error_Dec_2_Column == 'fixed':
        print('\nWARNING: Error_Dec_2 set to fixed value.')
        print('Fixed_Positional_Error = ', Fixed_Positional_Error)
        Error_Dec_2 = Fixed_Positional_Error
    else:    
        Error_Dec_2 = Data[Error_Dec_2_Column]

    print('\n')        
    Offset = 0.0 #IMPROVE THIS!!!
    Separation = Separation - Offset

    Positional_Error_Ra  = (Error_Ra_1**2  + Error_Ra_2**2) ** 0.5
    Positional_Error_Dec = (Error_Dec_1**2 + Error_Dec_2**2) ** 0.5
    Positional_Error = (Positional_Error_Ra + Positional_Error_Dec) / 2

    f = 1 / (2 * math.pi * Positional_Error**2) * math.e**(-Separation**2 / (2 * Positional_Error**2) )

    
    Width_Real = np.diff(Final_Real[:,0])[0]
    #print('Width_Real == ', Width_Real)
    #Real_Interpolation = np.interp(Magnitudes, Final_Real[:,0] + (Width_Real/2), Final_Real[:,1])
    Real_Interpolation_f = interp1d(Final_Real[:,0] + (Width_Real/2), Final_Real[:,1], kind='cubic')
    Real_Interpolation   = Real_Interpolation_f(Magnitudes)
    q_Normalised         = Real_Interpolation / np.sum(Final_Real[:,1]) * Q
    #Background_Interpolation = np.interp(Magnitudes, Final_Background[:,0] + (Width_Real/2), Final_Background[:,1])
    Background_Interpolation_f = interp1d(Final_Background[:,0] + (Width_Real/2), Final_Background[:,1], kind='cubic')
    Background_Interpolation   = Background_Interpolation_f(Magnitudes)
    n_Normalised               = Background_Interpolation / ( N_Original_First_Field * Match_Radius**2 * math.pi)

    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    plt.step(Final_Real[:,0], Final_Real[:,1], label='Final', where='post')
    plt.scatter(Magnitudes, Real_Interpolation, color='orange', s=0.2)
    plt.step(Final_Background[:,0], Final_Background[:,1], label='Final', where='post')
    plt.scatter(Magnitudes, Background_Interpolation, color='cyan', s=0.2)
    plt.xlabel('Magnitude')
    plt.ylabel('N')
    plt.axhline(y=0, color='black')
    plt.grid(True)
    plt.title('Interpolation')
    plt.savefig('log/PLOT_Interpolation.png')
    plt.close()



    Likelihood_Ratio = f * q_Normalised / n_Normalised
    Likelihood_Ratio_Column  = fits.Column(name='Likelihood_Ratio', array=Likelihood_Ratio, format='D')
    Likelihood_Ratio_Columns = fits.ColDefs([Likelihood_Ratio_Column])
    Data_Columns = Data.columns

    Merged_Columns = Data_Columns + Likelihood_Ratio_Columns

    Fits_File = fits.BinTableHDU.from_columns(Merged_Columns)
    Fits_File.writeto('Output/Matched_LR.fits', overwrite=True)
    print('Output saved into Output/Matched_LR.fits file.\n')

    Data_LR = fits.open('Output/Matched_LR.fits')[1].data
    print('Length of the Complete Matched Cat: ', len(Data_LR[ID_1]))
    
    Data_LR = Data_LR[ Data_LR[Magnitude_Column]>0 ]
    print('Value of Q used in LR calculations: ', Q)
    print('Length of the LR Column           : ', len(Data_LR[ID_1]))
    print('Length of the unique ID1 rows     : ', len(np.unique(Data_LR[ID_1]))  )
    print('Length of the LR > 0.2            : ', len(Data_LR [Data_LR['Likelihood_Ratio']>0.2 ][ID_1] ) )
    print('Length of the LR > 0.2 unique ID1 : ', len(np.unique(Data_LR [Data_LR['Likelihood_Ratio']>0.2 ][ID_1] )) )
   
    #Creating the log LR_Calculation output
    print('\nNormalisation Sum')
    print('np.sum(Final_Real[:,1]) == ', np.sum(Final_Real[:,1]), '\n')
    col1 = fits.Column(name='ID',                      format='20A', array=ID)
    col2 = fits.Column(name='Separation',              format='D',   array=Separation)
    col3 = fits.Column(name='Magnitudes',              format='D',   array=Magnitudes) 
    col4 = fits.Column(name='Real_Interpolation',      format='D',   array=Real_Interpolation)
    col5 = fits.Column(name='q_Normalised',            format='D',   array=q_Normalised)    
    col6 = fits.Column(name='Background_Interpolation',format='D',   array=Background_Interpolation)
    col7 = fits.Column(name='n_Normalised',            format='D',   array=n_Normalised)
    col8 = fits.Column(name='f',                       format='D',   array=f)
    col9 = fits.Column(name='Likelihood_Ratio',        format='D',   array=Likelihood_Ratio)
    cols= fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])    
    Fits_File = fits.BinTableHDU.from_columns(cols)
    Fits_File.writeto('log/LR_Calculation_Data.fits', overwrite=True)    

    
    return 








