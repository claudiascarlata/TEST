import numpy as np
cfrom astropy.io import ascii
#import pysynphot as S
import scipy as s
from scipy import optimize
from scipy.integrate import quad
import math as m
from scipy.interpolate import interp1d
import pyfits

def read_fits(name):
    datacube=pyfits.open(name)
    data=image[1].data
    err=image[2].data*5. +2
    L0=hdulist[0].header['LAMBDA0'] 
    DL=hdulist[0].header['DLAMBDA'] 
    N1=hdulist[0].header['NAXIS1'] 
    N2=hdulist[0].header['NAXIS2'] 
    out=['data':data,'err':err,'lambda0':L0,'deltal':DL,'n1':N1,'n2':N2]
    return out

read_fits('Par258_Obj146_G141.fits')
x=np.linspace(0,out['n2'])
for j in range(0,out['n1'])
    plt.plot(x,out['data'][:,j])

plt.show()


    
    



def compute_profile():
    load_SiII_1190()
    load_SiII_1193()
    load_SiII_1260()
    load_geom()
    velocity = np.linspace(-2500,2500,10000)
    dv=velocity[1]-velocity[0]
    #compute the convolution box to go to the correct spectral resolution
    box = int(30./dv)
    print 'Convolution box',box, dv
    c=3e5
    dx=c*(SiII_1193['resonant_wavelength']-SiII_1190['resonant_wavelength'])/SiII_1190['resonant_wavelength']
    shift = float(dx)
    idx_v = np.abs(velocity - dx).argmin()
    idx0 = np.abs(velocity - 0).argmin()
    dx =idx_v - idx0

    l_1190=compute_line(velocity,SiII_1190)
    l_1193=compute_line(velocity,SiII_1193)
    l_1260=compute_line(velocity,SiII_1260)
    em_l_1193=np.roll(l_1193['em_tot'],int(dx))
    abs_l_1193=np.roll(l_1193['abs'],int(dx))
#    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2,figsize=(10,8),tight_layout=True)
    fig, (ax1,ax2,ax3) = plt.subplots(3, 1,figsize=(7,15),tight_layout=True)

    ax1.plot(velocity,l_1190['tot'],'RoyalBlue',label='SiII 1190',linewidth=3)
    ax1.plot(velocity - shift,1.+em_l_1193-abs_l_1193+l_1190['fl_blu']+l_1190['fl_red'],'DarkOrange',label='SiII 1193',linewidth=3)
    ax1.plot(velocity,l_1260['tot'],'ForestGreen',label='SiII 1260',linewidth=3)
    ax1.legend(fontsize=15,loc=2)
    ax1.set_xlim(-499,200)
    ax1.set_ylim(0.0,2.2)
    ax1.set_xlabel('v [km/s]',fontsize=14)
    ax1.set_title('Original resolution')

    ax2.plot(velocity,1-l_1190['abs'],'RoyalBlue',linewidth=3)
    ax2.plot(velocity,1-l_1193['abs'],'DarkOrange',linewidth=3)
    ax2.plot(velocity,1-l_1260['abs'],'ForestGreen',linewidth=3)
    ax2.set_xlim(-499,200)
    ax2.set_ylim(0.0,1.2)
    ax2.set_xlabel('v [km/s]',fontsize=14)
    ax2.set_title('Only absorption')

    ax3.plot(velocity,smooth(l_1190['tot'],box)+ 0.2*np.random.randn(len(velocity)),'RoyalBlue',alpha=0.6)
    ax3.plot(velocity- shift,smooth(1.+em_l_1193-abs_l_1193+l_1190['fl_blu']+l_1190['fl_red'],box)+ 0.2*np.random.randn(len(velocity)),'DarkOrange',alpha=0.6)
    ax3.plot(velocity,smooth(l_1260['tot'],box)+ 0.2*np.random.randn(len(velocity)),'ForestGreen',alpha=0.6)
    ax3.set_xlim(-499,200)
    ax3.set_ylim(0.0,2)
    ax3.set_xlabel('v [km/s]',fontsize=14)    
    ax3.set_title('Simulated Observation')
 
    """    ax4.plot(velocity,smooth(l_1190['tot'],box)+ 0.2*np.random.randn(len(velocity)),'RoyalBlue',alpha=0.6)
    ax4.plot(velocity- shift,smooth(1.+em_l_1193-abs_l_1193+l_1190['fl_blu']+l_1190['fl_red'],box)+ 0.2*np.random.randn(len(velocity)),'DarkOrange',alpha=0.6)
    ax4.plot(velocity,smooth(l_1260['tot'],box)+ 0.2*np.random.randn(len(velocity)),'ForestGreen',alpha=0.6)
    ax4.set_xlim(-499,1500)
    ax4.set_ylim(0.0,2)
    ax4.set_xlabel('v [km/s]',fontsize=14)    
    ax4.set_title('Simulated Observation')
    """
    fig.savefig('discussion.png',dpi=400)
    fig.show()

    print load_geom()

