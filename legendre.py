#using legendre polynomial for fitting vega spectra data
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
import numpy.polynomial.legendre as nppl
#The data
data = np.loadtxt('vega.txt', unpack = True)

flux = data[1,:]
wave = data[0,:]
#Okay so basically we want to set up our rolling median and
#standard deviation here, it needs to roll to adjust to the
#craziness of our dataset
d = {'Wavelength' : wave, 'Intensity' : flux}

i_df = pd.DataFrame(d)
i_df['Std'] = i_df['Intensity'].rolling(50).std()
i_df['Median'] = i_df['Intensity'].rolling(50).median()
c_df = i_df.copy()

#Time for the function

def regression(df):
    x = np.arange(df.size)
    m,b,r,p,e = linregress(x,df)
    #liny = m * df['Wavelength'] + b
    yval = m * x[-1] + b
    return yval

#And it's time to initialize

a = 0.5
b = 1.0
iters = 6

for i in range(0,iters):
    upper = c_df['Median'] + a * c_df['Std']
    lower = c_df['Median'] - b * c_df['Std']

    mask = c_df['Intensity'] < upper
    mask = c_df['Intensity'] > lower

    c_df = c_df[mask]

    c_df['Std'] = c_df['Intensity'].rolling(50).std()
    c_df['Median'] = c_df['Intensity'].rolling(50).agg(regression)

#Now it's finally time for the Legendre fit baby

coefs = nppl.legfit(c_df['Wavelength'], c_df['Intensity'], 15)
fit = nppl.legval(c_df['Wavelength'], coefs)

plt.plot(wave, flux, 'bo')
plt.plot(i_df['Wavelength'], i_df['Intensity'])
plt.plot(c_df['Wavelength'], c_df['Intensity'])
plt.plot(c_df['Wavelength'], fit)
plt.title("Vega Spectra data with ${\sigma}$ clipped data")
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (W/m**2/m)')
plt.show()

plt.plot(wave,flux)
plt.title('Vega Spectra data')
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (W/m**2/m)')
plt.show()

# clipped = sigma_clip(flux, sigma_lower =2, sigma_upper = 4, iters=100)

# plt.plot(wave, clipped)
plt.show()
    
