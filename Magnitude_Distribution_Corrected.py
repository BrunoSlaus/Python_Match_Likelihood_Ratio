# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017

Function that, creates the magnitude distributon of the counterparts
for a matched field. It returns the data needed to plot the histograms
and the limit after which Corrected_Real becomes negative.
            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


def Magnitude_Distribution_Corrected(Matched_Field, Full_Counterpart_Field, Match_Radius, Magnitude_Column, N_Bins, Corrected_Background, Plot_Name, Range_Min, Range_Max ):
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import math

    #Histogram parameters
    Buffer    = 1      #Padding from min range to ensure "Limit" is calculated correctly


    
    Original_Field_Counts = len( np.unique(Matched_Field['Id']) )      #Number of matched sources in the orig. field
    #print('Original_Field_Counts: ', Original_Field_Counts)           #Printing the N of orig. matched sources
    
    Counterpart_Magnitudes = Matched_Field[Magnitude_Column]           #Getting magnitudes from the matched field
    Full_Magnitudes        = Full_Counterpart_Field[Magnitude_Column]  #Getting mag. from the full counterpart field

    Counterpart_Magnitudes = Counterpart_Magnitudes[np.logical_not(np.isnan(Counterpart_Magnitudes))] #Removing the NaN values
    Full_Magnitudes        =        Full_Magnitudes[np.logical_not(np.isnan(Full_Magnitudes))]        #Removing the NaN values

    #Creating the histograms using "np.histogram"
    Total_F, Total_Edges = np.histogram(Counterpart_Magnitudes, bins=N_Bins, range=(Range_Min,Range_Max))
    Width = np.diff(Total_Edges)
    Total_Edges = Total_Edges[:-1] #Removing the last unnecesarry edge
    Corrected_Background_F     = Corrected_Background[:,1]
    Corrected_Background_Edges = Corrected_Background[:,0]
    Corrected_Real_F           = Total_F - Corrected_Background_F
    Corrected_Real_Edges       = Total_Edges


    
    Total                = np.column_stack((  np.array(Total_Edges[:]), np.array(Total_F[:])  ))
    Corrected_Background = np.column_stack((  np.array(Corrected_Background_Edges[:]), np.array(Corrected_Background_F[:])  ))
    Corrected_Real       = np.column_stack((  np.array(Corrected_Real_Edges[:]), np.array(Corrected_Real_F[:])  ))



    #Finding the range where Corrected_Real is negative, and the limit of this regime.
    Corrected_Real_Minus       = Corrected_Real[Corrected_Real[:,1]<0]
    Corrected_Real_Minus_Edges = Corrected_Real_Minus[:,0]
    Corrected_Real_Minus_Edges = Corrected_Real_Minus_Edges[Corrected_Real_Minus_Edges>Range_Min+Buffer]  #Adding to eliminate the left min edge!!
    #print('q_Minus_Edges', Corrected_Real_Minus_Edges)
    if len(Corrected_Real_Minus_Edges)>0:
        Limit = np.amin(Corrected_Real_Minus_Edges)
    else:
        Limit = -99
        print('\nWARNING: No lower limit for negative Corrected_Real-value.')
    print('Calculatinge range where Corrected_Real<0. Printing the lower LIMIT:')
    print('Limit: ', Limit)
    print('\n')



    #Plotting the histograms
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    edges = Total_Edges           #All histogram edges are the same because the bins and min/max ranges do not differ

    plt.step(edges, Total_F, where='post', color="red")
    plt.step(edges, Corrected_Background_F, where='post', color="blue")
    plt.step(edges, Corrected_Real_F, where='post', color="orange")

    plt.xlabel('Magnitude')
    plt.ylabel('N')
    plt.axhline(y=0, color='black')
    plt.grid(True)
    plt.title('Magnitude Distribution ' + Plot_Name)
    #plt.vlines(Limit, -10, 30, color='red')
    plt.text(0.75, 0.75, ' RED=Real \n BLUE=Random \n ORANGE=Diff \n Match_Radius='+str(Match_Radius)+'\n Limit='+str(Limit), transform=ax.transAxes)
    plt.savefig('log/PLOT_Magnitude_Distribution_'+ Plot_Name +'.png')



    #Creating the function output
    col1 = fits.Column(name='Left-Bin_Edge', format='D', array=Total[:,0])
    col2 = fits.Column(name='Total',         format='D', array=Total[:,1])
    col3 = fits.Column(name='Background',    format='D', array=Corrected_Background[:,1])
    col4 = fits.Column(name='Real',          format='D', array=Corrected_Real[:,1])

    cols = fits.ColDefs([col1, col2, col3, col4])
    
    Fits_File = fits.BinTableHDU.from_columns(cols)
    Fits_File.writeto('log/Histogram_Data_Corrected.fits', overwrite=True)



    #Creating the function output 
    Function_Output = {'Total': Total, 'Corrected_Background': Corrected_Background, 'Corrected_Real': Corrected_Real, 'Limit': Limit}
    return Function_Output



