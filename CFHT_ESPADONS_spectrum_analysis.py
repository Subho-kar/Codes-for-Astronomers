#!/usr/bin/env python
#!/usr/anaconda3/envs/my_env_py3/bin python
import os
import os.path
import sys
import numpy as np
import warnings
import re
import shlex
import readline
import shutil
import subprocess
from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.constants import c
from astropy import units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
from astropy.convolution import Gaussian1DKernel, convolve
from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import equivalent_width, fwhm
from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)
import matplotlib.pyplot as plt

import pywt 

def LoadFitsSpectra(filename):
    
    hdul = fits.open(filename)
    print(hdul.info())
    Header=hdul[0].header
    Data=hdul[0].data
    wavedata=np.array(Data[0])
    fluxdata=np.array(Data[1])
    print("Zoom and decide the targeted line and close the window after that.")
    plt.figure()
    plt.plot(wavedata,fluxdata,'r')
    plt.ylabel('Relative flux')
    plt.xlabel('$\lambda (nm)$')
    plt.title('Whole spectrum')
    plt.grid() 
    plt.show()
    

    return wavedata,fluxdata


def AskAndLoadData(toask="Enter Filename : ", inpfilename=None):
    """ Asks user to enter filename and return the table loaded as numpy array. 
    toask is string which contains the prompt to ask."""
    
    while True:  #Keep asking till a file is properly loaded
        filename = input(toask).strip(' ') if inpfilename is None else inpfilename
        try :
            if os.path.splitext(filename)[1] == '.fits' : #User has inputed a fits file
                wavedata, fluxdata = LoadFitsSpectra(filename)
            else: 
                break
        except IOError :
            print("Error: Cannot find the file %s, Please give the correct file name."%(filename))
            inpfilename = None
        else:
            break


    return wavedata,fluxdata,filename

def Spec_extract(wavedata, fluxdata):
    #wavedata=np.array(Spec_data[0])
    #fluxdata=np.array(Spec_data[1])
    w1=0
    w2=0
    wavedata2=[]
    fluxdata2=[]
    w1, w2 = [float(x) for x in input("Enter the targeted wavelength range in integers(shorter first):").split()]
    indices=np.array(np.where((wavedata >= w1) & (wavedata < w2)))
    index=indices[:][0]
    """"Now we remove the extra rows by comparing the consective values and once the condition fails,
    we can break the loop"""
    for i in range(len(index)):
        
        
        if wavedata[index[i]] < wavedata[index[i+1]]:
           # print(index[i], wavedata[index[i]])
            wavedata2.append(wavedata[index[i]])
            fluxdata2.append(fluxdata[index[i]])
        else:
            break
    
    return wavedata2,fluxdata2

def Line_extract(toask="Enter line name:", inpline_name=None):
    """ Asks user to enter line name and return all the flux scaled data sets as a master astropy readable array. 
    toask is string which contains the prompt to ask."""
    
    PARENT_DIR=os.getcwd()  
    i=0
    ScaleSci=[]
    Scicorr=[]
    spec=[]
    Scidat_spec = []
    lam = []
    flam = []
    W = []
    FWHM = []
    while True:
        line_name = input("Enter the emission line (in nm) you want to analyse (format: HeI_447): ")
        file_path = os.path.join(PARENT_DIR, 'File_list_' + line_name)
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file1:
                filelist = [line.rstrip() for line in file1]
                for f in range(len(filelist)):
                    spec.append(np.loadtxt(filelist[f], delimiter=' '))         # Read from each file from the file list and store them as numpy arrays.
                Scicorr.append(spec[0][:,1])              
                for i in range(len(spec)-1):
                    ScaleSci.append(np.median(spec[0][:,1])/np.median(spec[i+1][:,1]))            
                    Scicorr.append(spec[i+1][:,1]*ScaleSci[i])
            break
        else:
            print(f"Error: Cannot find the file {file_path}. Please give the correct line name.")
            line_name = None
                                                                 
    lam.append(spec[0][:, 0] * u.nm)
    flam.append(Scicorr[0] * u.Unit('erg s-1 cm-2 AA-1'))
    Scidat_spec.append(Spectrum1D(spectral_axis=lam[0], flux=flam[0]))
    for i in range(1, len(spec)):
        lam.append(spec[i][:, 0] * u.nm)
        flam.append(Scicorr[i] * u.Unit('erg s-1 cm-2 AA-1'))
        Scidat_spec.append(Spectrum1D(spectral_axis=lam[i], flux=flam[i]))

    spec_region = SpectralRegion(lam[0][0], lam[0][-1])
    
    for j in range(len(Scidat_spec)):
        W.append(equivalent_width(Scidat_spec[j], regions=spec_region))
        FWHM.append(fwhm(Scidat_spec[j], regions=spec_region))
    
    return Scidat_spec, line_name, FWHM

def Stackspectra(Scidat_spec, FWHM):
    """ Loads the spectra into astropy readable table and computes the mean smoothed spectrum. Also returns  
        truncated data sets (in velocity scale) within the FWHM values."""
    
    Scidat_sum = Scidat_spec[0]
    for j in range(1, len(Scidat_spec)):
        Scidat_sum += Scidat_spec[j]
    
    for i in range(len(Scidat_spec)):
        plt.plot(Scidat_spec[0].spectral_axis, Scidat_spec[0].flux, label=f'Data {i+1}')
    plt.title('Scaled data sets at same flux level')
    plt.xlabel('$\lambda (nm)$')
    plt.ylabel('Relative flux')
    plt.legend()
    plt.grid()
    plt.show()
    Scidat_mean = Scidat_sum / len(Scidat_spec)
#    stddev = 60  # Standard deviation for Gaussian smoothing
    while True:
 #   Scidat_mean_gsmooth = convolve(Scidat_mean.flux, Gaussian1DKernel(stddev=stddev))
        val= int (input("Enter the value of standard deviation for the Gaussian function:"))
        Scidat_mean_gsmooth = gaussian_smooth(Scidat_mean, stddev=val)
        plt.plot(Scidat_mean.spectral_axis, Scidat_mean.flux, label='Mean spectrum')
        plt.plot(Scidat_mean_gsmooth.spectral_axis, Scidat_mean_gsmooth.flux, label='Mean Smoothed spectrum') 
        plt.xlabel('$\lambda (nm)$')
        plt.ylabel('Relative flux')
        plt.legend()
        plt.grid()
        plt.show()
        user_input = input('Are you satisfied with the smooth line profile? (yes/no): ')
        if user_input.lower() == 'no':
            continue
        elif user_input.lower() == 'yes':
            break
        else:
            print('Exiting because of wrong input!!!')  
            sys.exit()
    diff_spec = []
    diff_spec.append(Scidat_spec[0].flux - Scidat_mean_gsmooth.flux)
    
    for k in range(1, len(Scidat_spec)):
        diff_spec.append(Scidat_spec[k].flux - Scidat_mean_gsmooth.flux)
    
    """Select a particular wavelength range of the difference spectrum for WPS calculation."""
    w = float(input("Enter the central wavelength (in nm) of the spectral line being diagnosed: "))
    w1 = (w * u.nm) - FWHM[0] / 2
    w2 = (w * u.nm) + FWHM[0] / 2
    indices = np.where((Scidat_spec[0].spectral_axis >= w1) & (Scidat_spec[0].spectral_axis <= w2))
    indices = indices[0]
    red_spec = []
    table_1 = Table()
    print("In order to remove the edge effects, wavelength range to be considered between:", w1, " and ", w2) 
    for i in range(len(Scidat_spec)):
        table_2 = Table()
        wavedata3 = []
        fluxdata3 = []
        veldata = []
        for j in range(len(Scidat_spec[i].spectral_axis)):
            if (Scidat_spec[i].spectral_axis[j] >= w1) and (Scidat_spec[i].spectral_axis[j] <= w2):
                wavedata3.append(Scidat_spec[i].spectral_axis[j].value)
                fluxdata3.append((diff_spec[i][j].value))
                veldata.append((c * (wavedata3[-1] - w) / w).to(u.km/u.s).value)
        red_spec.append(np.column_stack((veldata, fluxdata3)))

    np.savetxt('Diff_spec_set-1', red_spec[0], delimiter = ' ')

    return diff_spec, red_spec
 

#def wavelet_filter(red_spec, line_name):
#    """ Filter the spectra based upon its various scaling components using a mother wavelet function """

#    while True:
#        filtspec = []
#        recon = []
#        print("The available mother wavelet functions are:", pywt.wavelist(family=None, kind='discrete')) 
#        wavelet_type, a = input("Enter the desired type of mother wavelet function (eg: db6) and decomposition level number (any integer > 0) :").split()
#        lvl = int(a) 
#        for i in range(len(red_spec)):
#            coeff=pywt.wavedec(red_spec[i][:,1],wavelet_type,mode='symmetric', level=lvl, axis=-1)        #Multi-Wavelet filter
#            for j in range(5, len(coeff)):
#                coeff[j] = np.zeros_like(coeff[j])
#            #coeff[-1:-5] = np.zeros_like(coeff[-1:-5])                                                   #Ignore the last two coefficient array
#            recon=pywt.waverec(coeff,wavelet_type,mode='symmetric', axis=-1)                #Multilevel 1D Inverse Discrete Wavelet Transform
#            filtspec.append(recon)
            
#        for i in range (len(filtspec)):
#            plt.plot(red_spec[0][:,0],filtspec[i][0:len(filtspec[i])-1] + i , label=f'Data {i+1}')
#        plt.title(f'{line_name} Filtered spectra')
#        plt.xlabel('$velocity (km/s)$')
#        plt.ylabel('Data sets (at 40s interval)')
#        plt.grid()  
#        plt.show() 
#        user_input = input('Are you satisfied with the wavelet transformed spectra? yes/no: ')
#        if user_input.lower() == 'no':
#            continue
#        elif user_input.lower() == 'yes':
#            break
#        else:
#            print('Exiting because of wrong input!!!')  
#            sys.exit()
#    return filtspec   
 

#def wave_pow_spec(filtspec):
#    """ Wavelet power spectrum calculation """
#    power_spectrum = []
#    tot_pow_spec = np.zeros_like(filtspec[0])  # Initialize tot_pow_spec as zeros array
#    avg_pow_spec = np.zeros_like(filtspec[0])  # Initialize avg_pow_spec as zeros array
    
#    for i in range(len(filtspec)):
#        power_spectrum.append(np.abs(filtspec[i]) ** 2)
#        tot_pow_spec += power_spectrum[i]
    
#    avg_pow_spec = tot_pow_spec / len(filtspec)

#    return avg_pow_spec


def main():
  
    PARENT_DIR=os.getcwd()
    line_name = None
    print('Working on the directory : ', PARENT_DIR)
    while True:
        
        #directories=open(os.path.join(,PC.OUTDIR,'directories'),'w')
        scifname = None     #File initialisation
        wavedata, fluxdata, scifname = AskAndLoadData(toask="Enter the name of science star spectrum file :", inpfilename=scifname)
        wavedata2, fluxdata2 = Spec_extract(wavedata, fluxdata)
        plt.plot(wavedata2,fluxdata2)
        plt.xlabel('$\lambda (nm)$')
        plt.ylabel('Relative flux')
        plt.grid()
        plt.show()
        Specline_data=np.column_stack((wavedata2, fluxdata2))
       #For extracting the output spectrum in ascii format 
        fnameSuffix=input("Enter the name of the emission line that was extracted(format: HeI_447): ")
        scifnameNew_ascii = os.path.splitext(scifname)[0]+'_'+fnameSuffix
        file_path=os.path.join(PARENT_DIR, fnameSuffix)
        print(file_path)  

        if os.path.exists(file_path):                               
            os.chdir(file_path)
        else:                                                      
            os.mkdir(file_path)
            os.chdir(file_path)  
        
        np.savetxt(scifnameNew_ascii,Specline_data, delimiter=' ')
        print('The spectra in ascii format was saved by the following name:',scifnameNew_ascii) 
        with open(os.path.join(PARENT_DIR,'File_list_'+fnameSuffix),'a') as file1:         #Save the ascii filenames in a filelist
            file1.write(scifnameNew_ascii + '\n')
            print('File name saved in', os.path.join(PARENT_DIR,'File_list_'+fnameSuffix))          
        user_input = input('Do you want to continue the same for other spectra? yes/no: ')
        if user_input.lower() == 'yes':
            os.chdir(PARENT_DIR)
            continue
        elif user_input.lower() == 'no':
            break
        else:
            print('Exiting because of wrong input!!!')  
            break

#    Scidat_spec, line_name, FWHM = Line_extract(toask="Enter the emission line (in nm) you want to analyse (format: HeI_447):", inpline_name=line_name)
    #print(FWHM)
    #diff_spec, red_spec =  Stackspectra(Scidat_spec, FWHM)
#    diff_spec, red_spec =  Stackspectra(Scidat_spec, FWHM)
#    print("Close the window before proceeding with wavelet filtering")
#    print("!!!Save the plots before closing the windows!!!")
#    for i in range (len(diff_spec)):
#        plt.plot(Scidat_spec[0].spectral_axis, diff_spec[i] + i * u.Unit('erg s-1 cm-2 AA-1') , label=f'Set {i+1}')
#    plt.title(f'{line_name} Difference spectra')
#    plt.xlabel('$\lambda (nm)$')
#    plt.ylabel('Data sets (at 40s interval)')
#    plt.legend()
#    plt.grid()  
#    plt.show()
#    for i in range (len(red_spec)):
#        plt.plot(red_spec[i][:,0], red_spec[i][:,1] + i,'.', label=f'Set {i+1}')
#    plt.title(f'Windowed difference spectra of {line_name}')
#    plt.xlabel('$Velocity (km/s)$')
#    plt.ylabel('Data sets (at 40s interval)')
#    plt.legend()
#    plt.grid()  
#    plt.show()
#    filt_spec = wavelet_filter(red_spec, line_name)
#    avg_pow_spec = wave_pow_spec(filt_spec)
#    #print(avg_pow_spec.ndim)
#        raise ValueError("avg_pow_spec should be a 2D array")
#    plt.plot(red_spec[0][:,0], np.log10(avg_pow_spec[0:len(avg_pow_spec)-1]), '-r')
#    plt.title("Average Wavelet Power Spectrum")
#    plt.xlabel("Scale (km/s)")
#    plt.ylabel("log(Power)")
#    plt.grid()
#    plt.show()
#    print('***********Analysis was completed************')
if __name__ == "__main__":
    main() 
