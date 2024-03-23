import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, widgets
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

## Here, we are defining some parameters
G = 6.67e-11 #units are m^3 kg^-1 s^-2
M_sun = 1.99e30 #units is kg
au = 1.5e11 # the distance of the earth to the sun in m
ly = 9.46e15 # lightyear in m
c = 2.998e8 # speed of light in m/s 

def kepler_velocity(a, M):
    """
    This function calculated the velocity given a semi-major axis a and mass of the central object M.
    a: semi-major axis in meters
    M: central mass in kg
    returns: velocity in m/s
    """
    v = np.sqrt(G*M/a)
    return v

def schwarzschild(Mbh):
    """
    This function calculated the Schwarzschild radius of a black hole
    MBh: black hole mass in kg
    returns: Schwarzschild radius in m
    """
    rs = 2*G*Mbh/(c**2)
    return rs

def bh_interactive(manual=True, log=False):
    #setting the scale over which we will look at the rotation velocity.
    #min_r = 0.01
    if log:
        max_r = 1000
    else:
        max_r = 100

    ##We will not be looking at distances in meters, instead, lets use astronimcal units
    scale_by = au
    scale_by_name = 'Astronomical Units'
    

    def f(r_blr, bhmass):
        #creating the figure
        plt.xkcd()
        plt.figure(2)
        ##Calculate the schwarzschild radius
        r_s = schwarzschild(10**(bhmass)*M_sun)/scale_by
        ##creates an array of distances
        if log:
            #x = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,50,100,500,1000])
            x = np.logspace(np.log10(r_s), np.log10(max_r), num=1000)
        else:
            x = np.linspace(r_s, max_r, num=1000)
        ## calculates keplerian velocity
        v_ar = kepler_velocity(x*scale_by, 10**(bhmass)*M_sun)/1000
        ##plots the Keplerian rotation curve
        plt.plot(x[x<3*r_s], v_ar[x<3*r_s], c='grey', ls='--')
        plt.plot(x[x>3*r_s], v_ar[x>3*r_s], c='k')
        ##calculates and markes the current radius and velocity
        plt.axvline(r_blr, c='m', ls='--')
        if r_blr >= r_s:
            v = kepler_velocity(r_blr*scale_by, 10**(bhmass)*M_sun)/1000
        else:
            v = kepler_velocity(r_s*scale_by, 10**(bhmass)*M_sun)/1000
        plt.axhline(v, c='m', ls='--')
        plt.plot([r_blr], [v], marker='*', c='m', ms=20)
        ##Calculate the schwarzschild radius
        r_s = schwarzschild(10**(bhmass)*M_sun)/scale_by
        #Mark the innermost stable circular orbit
        plt.text(3.1*r_s, 0.9*max(v_ar), "No more stable circular orbits", rotation=90, ha='left', va='top', c='orange', weight='bold')
        plt.axvline(3*r_s, c='orange')
        plt.axvline(r_s, c='orange')
        if log is False:
            plt.gca().add_patch(Rectangle((r_s,0),2*r_s,1.2*max(v_ar),fill=True, color='orange', alpha=0.5, zorder=100))
        #Mark the schwarzschild radius
        plt.text(1.1*r_s, 0.9*max(v_ar), "Schwarzschild Radius", rotation=90, ha='left', va='top', c='r', weight='bold')
        plt.axvline(r_s, c='r')
        plt.axvline(0, c='r')
        if log is False:
            plt.gca().add_patch(Rectangle((0,0),r_s,1.2*max(v_ar),fill=True, color='r', alpha=0.5, zorder=100))
        if r_blr > 3*r_s:
            plt.text(0.25*max_r, 0.8*max(v_ar),
                     "You are at %.2f %s \n from a %.2f Million Solar mass black hole \n and are rotating at %.2f km/s"
                     %(r_blr, scale_by_name, 10**(bhmass)/1e6, v), weight='bold')
        elif r_blr > r_s:
            plt.text(0.25*max_r, 0.8*max(v_ar),
                     "You are no longer on a stable orbit and \n are plunging towards the black hole!", weight='bold', c='orange')
        else:
            plt.text(0.25*max_r, 0.8*max(v_ar),"You have fallen into the black hole!", weight='bold', c='r')
        if log:
            plt.xscale("log")
            plt.xlim(0.5*r_s, max_r)
        else:
            plt.xlim(0,max_r)
        plt.xlabel("Distance from black hole (in %s)" % scale_by_name)
        plt.ylabel("Orbital Velocity (in km/s)")
        plt.title("Orbiting a black hole")
        plt.show()
    
    
    if log:
        blr_slider = widgets.FloatText(
                    value=7.5,
                    description='distance',
                    disabled=False)
    else:
        blr_slider = widgets.FloatSlider(
             value=0.5*max_r,
             min=0,
             max=max_r,
             step=0.1,
             description='distance',
             disabled=False,
             continuous_update=True,
             orientation='horizontal',
             readout=True,
             readout_format='.1f')

    bhmass_slider = widgets.FloatSlider(
        value=7,
        min=6,
        max=9,
        step=0.1,
        description='(log$M_{BH}$)',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')

    interactive_plot = interactive(f, {'manual': manual}, r_blr=blr_slider, bhmass=bhmass_slider)
    return(interactive_plot)


def analyze_spectrum(manual=True):
    wl = np.arange(6950, 7400, 5)
    def emlinecont(cont, linepos, velocity, lineflux):
        sigma = (velocity*linepos)/(c/1000)
        norm = lineflux*1000
        line = (norm / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((wl - linepos)** 2)/ (2 * sigma**2))
        full = line+cont
        return(full)
    spectrum = np.loadtxt('spectrum.txt', delimiter=',')
    wavelength = spectrum[:,0]
    flux = spectrum[:,1]
    
    def f(cont, linepos, velocity, lineflux):
        #creating the figure
        plt.xkcd()
        plt.figure(2)
        plt.plot(wavelength, flux, ls='-', c='k', lw=1)
        plt.xlabel('Wavelength [$\AA$]')
        plt.ylabel('Flux')
        plt.xlim(6950, 7400)
        plt.ylim(12, 27)
        z=0.4563313
        plt.title("Measuring gas velocity around a black hole.")
        plt.axvline(5008*(1+z), ls=':', c='b')
        plt.axvline(4960*(1+z), ls=':', c='b')
        plt.axvline(4863*(1+z), ls='--', c='r')
        plt.plot(wl, emlinecont(cont, linepos, velocity, lineflux), ls='-', c='r')
        
    cont_slider = widgets.FloatSlider(
        value=15,
        min=10,
        max=20,
        step=1,
        description='baseline',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')

    linepos_slider = widgets.FloatSlider(
        value=7100,
        min=6950,
        max=7400,
        step=1,
        description='position',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')
    
    velocity_slider = widgets.FloatSlider(
        value=1000,
        min=100,
        max=10000,
        step=10,
        description='velocity',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')
    
    lineflux_slider = widgets.FloatSlider(
        value=1,
        min=0.1,
        max=5,
        step=0.1,
        description='lineflux',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')


    interactive_plot = interactive(f, {'manual': manual}, cont=cont_slider, linepos=linepos_slider,
                                  velocity=velocity_slider, lineflux=lineflux_slider)
    return(interactive_plot)



def fit_spectrum(manual=True):
    def emlinecont(cont, linepos, velocity, lineflux):
        sigma = (velocity*linepos)/(c/1000)
        norm = lineflux*1000
        line = (norm / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((wl - linepos)** 2)/ (2 * sigma**2))
        full = line+cont
        return(full)
    spectrum = np.loadtxt('spectrum.txt', delimiter=',')
    full_wavelength = spectrum[:,0]
    full_flux = spectrum[:,1]
    wlmask = (full_wavelength > 6950) & (full_wavelength < 7200)
    wl = full_wavelength[wlmask]
    flux = full_flux[wlmask]
    z=0.4563313
        
    def f(cont, linepos, velocity, lineflux):
        #creating the figure
        linefit = emlinecont(cont, linepos, velocity, lineflux)
        plt.xkcd()
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(20,10))
        ax0.plot(wl, flux, ls='-', c='k', lw=1)
        ax0.set_ylabel('Flux')
        ax0.set_xlabel("Wavelength [\AA]")
        ax0.set_ylim(12, 18)
        ax0.set_title("Spectrum and Fit")
        z=0.4563313
        ax0.axvline(4863*(1+z), ls='--', c='r')
        ax0.plot(wl, linefit, ls='-', c='r')
        
        ax1.set_title("Distance beween spectrum and fit")
        ax1.set_xlabel("Wavelength [\AA]")
        ax1.plot(wl, flux-linefit, ls='-', c='k')
        ax1.axhline(0,ls='--', c='grey')
        
    cont_slider = widgets.FloatSlider(
        value=12,
        min=10,
        max=20,
        step=1,
        description='baseline',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')

    linepos_slider = widgets.FloatSlider(
        value=7100,
        min=6950,
        max=7400,
        step=1,
        description='position',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')
    
    velocity_slider = widgets.FloatSlider(
        value=1000,
        min=100,
        max=10000,
        step=10,
        description='velocity',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f')
    
    lineflux_slider = widgets.FloatSlider(
        value=1,
        min=0.01,
        max=2,
        step=0.01,
        description='lineflux',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.2f')


    interactive_plot = interactive(f, {'manual': manual}, cont=cont_slider, linepos=linepos_slider,
                                  velocity=velocity_slider, lineflux=lineflux_slider)
    return(interactive_plot)

def tellmetheblackholemass():
    print("The mass of the black hole is 8 x 10^8 solar masses or 1.6 x 10^39 kg.")