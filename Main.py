# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017


            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
from Topcat_Match import skymatch
from Topcat_Match_Tight import skymatch_tight
from Topcat_Subset_Inf import Topcat_Subset_Inf
from Distribution_Merger import Distribution_Merger
from Magnitude_Distribution import Magnitude_Distribution
from Magnitude_Distribution_Corrected import Magnitude_Distribution_Corrected
from Random_Sample import Random_Sample
from Background_Corrected import Background_Corrected
from Likelihood_Ratio import Likelihood_Ratio
from Final_Plotter import Final_Plotter
from Reliability import Reliability
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import astropy.units as u
from subprocess import call


############################################################################
#                    Sliding Parameters:                                   # 
############################################################################
Input_Folder     = 'Input/'      #Name of the folder with the input fields
Output_Folder    = 'Output/'     #Name of the folder where the matched field is stored

Field_Name_1     = 'XXL-N_gmrt_out_CorrIRACch1.fits'
Field_Name_2     = 'XXL-N_irac_ch1.fits'

Match_Radius     = 3             #Match radius in arcsec

Magnitude_Column   = 'MAG_AUTO_ch1_IRAC'
ID_1               = 'ID'
ID_2               = 'NUMBER'
Error_Ra_1_Column  = 'E_Ra'
Error_Dec_1_Column = 'E_Dec'
Error_Ra_2_Column  = 'fixed'
Error_Dec_2_Column = 'fixed'

Full_Area        = 16.7*3600*3600   #Area of the full counterpart field in arcsec cca 60??
Force_Q          = 'No'
Q_Forced         = 0.8 

N_Bins    = 40
Range_Min = 8    #Histogram min range
Range_Max = 28    #Histogram max range

Move_Rescale = 0
Match_Rescale= 12/3
Sample_Match_Divide = 6 #Dividing the radius for sample mag distribution match

Force_Merge_Limit        ='No' #For merging the histogrms
Forced_Merge_Limit_Value = 21   #Value

N_Sample = 5000      #Sample lenth

Set_N_Original_First_Field   = 'Yes'
N_Original_First_Field_Value = 4615
############################################################################
print('\n******************************')
print('Starting the Main.py code.')
print('******************************\n')

Field_1 = Input_Folder+Field_Name_1
Field_2 = Input_Folder+Field_Name_2
print('Field_1      = ',Field_1)
print('Field_2      = ',Field_2)
print('Match_Radius = ',Match_Radius)
print('Move_Rescale = ',Move_Rescale)
print('Match_Rescale= ',Match_Rescale)

#We need only the part of the counterpart field with the nonzero Magnitude_Column
Full_Counterpart_Field = fits.open(Field_2)[1].data
Full_Counterpart_Field_NoNan  = Full_Counterpart_Field[np.logical_not(np.isnan(Full_Counterpart_Field[Magnitude_Column]))]
Full_Counterpart_Field_NoNeg  = Full_Counterpart_Field_NoNan[ Full_Counterpart_Field_NoNan[Magnitude_Column]>0 ]
Full_Counterpart_Field_Length = len(Full_Counterpart_Field_NoNeg[ID_2])
print('Length of the Full Magnitude Column : ', Full_Counterpart_Field_Length)

Topcat_Subset_Inf(Field_2, Magnitude_Column, 0, 'Infinity', 'log/Full_Counterpart_Field_NoNegMag.fits')
Field_Name_2 = 'Bla'
Field_2      = 'log/Full_Counterpart_Field_NoNegMag.fits'



if Set_N_Original_First_Field == 'Yes':
    N_Original_First_Field =N_Original_First_Field_Value
    print('\nWARNING: Setting N_Original_First_Field manually.')
else:    
    Original_First_Field   = fits.open(Field_1)[1].data
    N_Original_First_Field = len(  np.unique(Original_First_Field['Id'])  ) 
print('N_Original_First_Field= ',N_Original_First_Field)

#Matching the two fields
print('\n\n********************\nPART 1: MATCHING THE FIELDS:')
skymatch(
    Field_1,
    ["Ra", "Dec"],
    Field_2,
    ["Ra", "Dec"],
    Match_Radius, Output_Folder+'Match.fits')



#Measuring the initial magnitude distribution
print('\n\n********************\nPART 2: INITIAL MAGNITUDE DISTRIBUTION:')
Matched_Field   = fits.open(Output_Folder+'Match.fits')[1].data
N_Matched_Field = len(Matched_Field[Magnitude_Column])

Full_Counterpart_Field = fits.open(Field_2)[1].data
print('Length of the Full Counterpart Field (with Maggnitude_Column) : ', len(Full_Counterpart_Field[Magnitude_Column]) )

Initial_Distribution = Magnitude_Distribution(Matched_Field, Full_Counterpart_Field, Match_Radius, Magnitude_Column, N_Bins, Full_Area, 'Initial', Range_Min, Range_Max, N_Original_First_Field )    
Initial_Total        = Initial_Distribution['Total']
Initial_Real         = Initial_Distribution['Real']
Initial_Background   = Initial_Distribution['Background']
Limit                = Initial_Distribution['Limit']

print('\nCalculating the Q estimation:')
print('N_Original_First_Field == ', N_Original_First_Field)
print('N_Matched_Field        == ', N_Matched_Field)
print('Initial_Total SUM      == ', np.sum(Initial_Total[:,1]) )
print('Match_Radius           == ', Match_Radius)
print('Initial_Background SUM == ', np.sum(Initial_Background[:,1]) )
Q_Up    = N_Matched_Field - (np.sum(Initial_Background[:,1]) )
Q       = Q_Up / N_Original_First_Field
print('Value of Q parameter   == ', Q)
print('Initial_Real SUM       == ', np.sum(Initial_Real[:,1]) )
Q_Check = np.sum(Initial_Real[:,1]) / N_Original_First_Field
print('Checking Q parameter   == ', Q_Check)
if Force_Q == 'Yes':
    Q = Q_Forced
    print('Forcing Q to value ', Q)
else:
    print('Using the calculated Q value of ', Q)


#MATCHING FOR SAMPLE i.e. TIGHT MATCH
print('\n\n********************\nPART 3: TIGHT MATCH:')
skymatch_tight(
    Field_1,
    ["Ra", "Dec"],
    Field_2,
    ["Ra", "Dec"],
    Match_Radius/Sample_Match_Divide, Output_Folder+'Tight_Match.fits')

#Measuring the SAMPLE magnitude distribution
Tight_Matched_Field = fits.open(Output_Folder+'Tight_Match.fits')[1].data
print('Length of the Random-sample         : ', N_Sample)

Tight_Matched_Magnitudes = Tight_Matched_Field[Magnitude_Column]
Tight_Matched_Magnitudes = Tight_Matched_Magnitudes[np.logical_not(np.isnan(Tight_Matched_Magnitudes))]

Tight_Matched_F, Tight_Matched_Edges = np.histogram(Tight_Matched_Magnitudes, bins=N_Bins, range=(Range_Min,Range_Max))
Tight_Matched                        = np.column_stack((  np.array(Tight_Matched_Edges[:-1]), np.array(Tight_Matched_F[:])  ))



#Creating a random sample of sources from the full counterpart field
print('\n\n********************\nPART 4: CREATING A SUBSAMPLE:')
Sample = Random_Sample(Field_2, N_Sample, Tight_Matched, Magnitude_Column, Match_Radius, Move_Rescale)
Random_Sample_Data = Sample['Sample_Array']
Random_Sample_Name = Sample['Random_Sample_Name']



#Measuring the magnitude distribution of THE SAMPLE
print('\n\n********************\nPART 5: CORRECTED BACKGROUND CALCULATION:')
Corrected_Background   = Background_Corrected(Random_Sample_Name, Field_2, Match_Radius, Initial_Background, Magnitude_Column, N_Bins, Range_Min, Range_Max, N_Sample, N_Original_First_Field, Match_Rescale)['Corrected_Background']
Corrected_Distribution = Magnitude_Distribution_Corrected(Matched_Field, Full_Counterpart_Field, Match_Radius, Magnitude_Column, N_Bins, Corrected_Background, 'Corrected', Range_Min, Range_Max )    
Corrected_Real         = Corrected_Distribution['Corrected_Real']



#Obtaining the final q distribution
print('\n\n********************\nPART 6: CREATING THE FINAL MAGNITUDE DISTRIBUTION:')
Final_Real_Distribution = Distribution_Merger(Initial_Real, Corrected_Real, Limit, 'Real', Force_Merge_Limit, Forced_Merge_Limit_Value)['Final_Distribution']
Final_Background        = Distribution_Merger(Initial_Background, Corrected_Background, Limit, 'Background', Force_Merge_Limit, Forced_Merge_Limit_Value)
Final_Background_Distribution = Final_Background['Final_Distribution']
Final_Background_Limit        = Final_Background['New_Limit']
print('Final_Background_Limit ', Final_Background_Limit)



#Plotting the final distribution
print('\n\n********************\nPART 7: PLOTTING THE FINAL MAGNITUDE DISTRIBUTION:')
Final_Plotter(Matched_Field, Final_Real_Distribution, Final_Background_Distribution, Final_Background_Limit, Match_Radius, Magnitude_Column, N_Bins, 'FINAL', Range_Min, Range_Max)



#Calculating the LR & REL of each source
print('\n\n********************\nPART 8: CREATING THE LR & REL COLUMNS:')
Likelihood_Ratio(Output_Folder+'Match.fits', Magnitude_Column, Final_Real_Distribution, Final_Background_Distribution, ID_1, Error_Ra_1_Column, Error_Dec_1_Column, Error_Ra_2_Column, Error_Dec_2_Column, N_Original_First_Field, Match_Radius, Q)
Reliability(Q)



call('rm -f log/Full_Counterpart_Field_NoNegMag.fits', shell=True)  #REMOVING THE Field 
print('\n******************************')
print('FINISH: Ending the code.')
print('******************************\n')

