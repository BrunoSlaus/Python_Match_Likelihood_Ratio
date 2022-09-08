# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Fri May 12 16:03:09 2017

Function that..
            
Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


def Distribution_Merger_STD(Distribution_1, Distribution_2, Distribution_2_STD, Limit, Plot_Name, Force_Merge_Limit, Forced_Merge_Limit_Value):
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import math    

    if Force_Merge_Limit=='Yes':
        Limit = Forced_Merge_Limit_Value
        print('\nWARNING:Forcing Limit value during distribution merging.')
        print('Forced_Merge_Limit_Value = ', Forced_Merge_Limit_Value, '\n')



  
    #Opening the data
    Edges    = Distribution_1[:,0]
    Width    = np.diff(Edges)
    Width    = np.append(Width, Width[0])
    Values_1 = Distribution_1[:,1]
    Values_2 = Distribution_2[:,1]
    Values_2_STD = Distribution_2_STD[:,1]
    #print(Edges)
    #print(Values_1)
    #print(len(Edges))
    #print(len(Values_1))
    
    #Plotting the data
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    plt.step(Edges+Width, Values_1, label='Initial', where='pre')
    plt.step(Edges+Width, Values_2, label='Corrected', where='pre')
    plt.errorbar(Edges+(Width/2), Values_2, yerr=Values_2_STD, xerr=None, fmt='o', markersize=0.5)
    plt.xlabel('Magnitude')
    plt.ylabel('N')
    plt.title(Plot_Name +' Distribution Merged')
    #plt.vlines(Limit, -10, 30, color='red')
    #plt.text(0.75, 0.75, ' RED=Initial \n BLUE=Initial \n Limit='+str(Limit), transform=ax.transAxes)
    plt.savefig('log/'+ Plot_Name +'_Distribution_Merged.png')
    plt.close()

    #Merging the distributions
    Left_Part = Distribution_1[Distribution_1[:,0]<Limit]
    Right_Part= Distribution_2[Distribution_2[:,0]>=Limit]
    Final_Distribution = np.concatenate([Left_Part, Right_Part])

    #Plotting the final distribution
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    plt.step(Edges+Width, Final_Distribution[:,1], label='Final', where='pre')
    plt.xlabel('Magnitude')
    plt.ylabel('N')
    plt.title(Plot_Name +' Distribution Final')
    #plt.vlines(Limit, -10, 30, color='red')
    plt.savefig('log/'+ Plot_Name +'_Distribution_Final.png')
    plt.close()


    return {'Final_Distribution':Final_Distribution, 'New_Limit': Limit}    #The function returns the radius and the number of matches
