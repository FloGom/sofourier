# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 2021

@author: flogom
"""

__version__='0.1'

# === Imports === #
# 1. std libraries import
from os import mkdir, path, environ

# 2. third party libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

# 3. user defined libraries, functions

def len_2pts(pt1, pt2):
    res = np.sqrt((pt2[1] - pt1[1])**2 + (pt2[0] - pt1[0])**2)
    return(res)

def total_lenght(pts):
    nb_pts = len(pts)
    res = 0
    for i in range(1, nb_pts):
        res += len_2pts(pts[i], pts[i-1])
    return(res)

def parametric_lin(pts):

    pts_array = np.array(pts)

    x_val = pts_array[:,0]
    y_val = pts_array[:,1]

    nb_points = len(pts)
    # tot_len = total_lenght(pts)
    t = [0]

    for i in range(1, nb_points):
        my_len = len_2pts(pts[i-1], pts[i])
        delta_t = my_len
        t += [t[i-1] + delta_t]
    
    return(t, x_val, y_val)

def four_cn(t, y_t, n):
    "compute the complex fourier coefficient at rank n"
    period = max(t)
    if n==0:
        cout = y_t.mean()
    else:
        c = y_t*np.exp(-1j*2*n*np.pi*t/period)
        cout = c.sum()/c.size
    return(cout)

def four_ccoeff(t, y_t, max_n):
    "compute the n first complex fourier coefficient"
    coeffs = [four_cn(t, y_t, n) for n in range(max_n+1)]
    return(coeffs)

def ifour_ccoeff(t, coeffs, max_n):
    "compute the inverse fourier transform with the first max_n ranks complex coeffs"
    period = max(t)
    y_t_if = []
    for cur_t in t:
        y = np.array([2*coeffs[i]*np.exp(1j*2*i*np.pi*cur_t/period) for i in range(1,max_n+1)])
        y_t_if += [y.sum()]
    y_t_if = np.array(y_t_if)
    y_t_if += coeffs[0]
    return(y_t_if)

class Sofourier(object):
    ""
    def __init__(self, file_path, **kwargs):
        """
        file_path : string
            path to the .csv file holding data points X, Y
        kwargs : dict
            kwarg from pd.read_csv
        """
        self.fpath = file_path
        with open(file_path, 'rb') as myf:
            df_data = pd.read_csv(myf, **kwargs, header=None)
        
        x_val = df_data.iloc[:,0]
        y_val = df_data.iloc[:,1]
        nb_pts = len(x_val)
        xy_closedpath = []
        xy_closedpath += [[x_val[i], y_val[i]] for i in range(nb_pts)] 
        xy_closedpath += [xy_closedpath[0]]
        xy_closedpath = np.array(xy_closedpath)
        self.x_val = xy_closedpath[:,0]
        self.y_val = xy_closedpath[:,1]
        self.nb_pts = len(x_val)
        #parametrization
        t, x_t, y_t = parametric_lin(xy_closedpath)
        self.t_rspl = np.linspace(0, max(t), 1000)
        self.x_t_rspl = np.interp(self.t_rspl, t, x_t)
        self.y_t_rspl = np.interp(self.t_rspl, t, y_t)

        # center of mass
        self.x_mean = self.x_t_rspl.mean()
        self.y_mean = self.y_t_rspl.mean()

        self.bcolor = 'black'
        self.mcolorh = 0.61
        self.cmode = 'maincolor'
        self.dcolor = 'white'

        self.maxrank = 15

        self.min_sat = 0
        self.max_sat = 1

        self.min_alpha = 0.2
        self.max_alpha = 1


    def set_color(self, color_mode='maincolor', main_color_hue=None, data_color='white', background='black'):
        """
        color_mode : str
            'auto' : rainbow style
            'maincolor' : one color defined by its hue (0 to 1)
        main_color_hue : float
            hue code number between 0 and 1
        data_color : matplotlib color
            color to use for the data curve, usually white or black
        background : matplotlib color
            background color of the fig
        """
        self.bcolor = background
        self.mcolorh = main_color_hue
        self.cmode = color_mode
        self.dcolor = data_color
    
    def set_maxrank(self, max_rank):
        "Set the maximum number of Fourier rank to draw, starting from 1"
        self.maxrank = max_rank
    
    def set_sat(self, min_sat, max_sat):
        """
        Set saturation levels for the 'maincolor' mode, starting from center of
        figure to outer curve.
        """
        self.min_sat = min_sat
        self.max_sat = max_sat

    def set_alpha(self, min_alpha, max_alpha):
        """
        Set alpha levels (all color modes), starting from outer curve to center 
        curve
         """
        self.min_alpha = min_alpha
        self.max_alpha = max_alpha

    def plot(self):
        "Plot the Sofourier fig"
        # prepare data
        x_ccoeff = four_ccoeff(self.t_rspl, self.x_t_rspl, self.maxrank)
        y_ccoeff = four_ccoeff(self.t_rspl, self.y_t_rspl, self.maxrank)

        # plot data
        fig, ax = plt.subplots()
        fig.set_facecolor(self.bcolor)
        ax.set_aspect('equal')
        plt.axis('off')
    
        delta_s = (self.max_sat - self.min_sat)/(self.maxrank-1)
        delta_a = (self.max_alpha - self.min_alpha)/(self.maxrank-1)

        for r in range(1, self.maxrank):
            if self.cmode == 'auto':
                mycolor = hsv_to_rgb(((r-1)/(self.maxrank-1), self.max_sat-(r-1)*delta_s, 1))
            else:
                mycolor = hsv_to_rgb((self.mcolorh, self.max_sat-(r-1)*delta_s, 1))
            
            myalpha = self.min_alpha + (r-1)*delta_a

            x_if = ifour_ccoeff(self.t_rspl, x_ccoeff, r)
            y_if = ifour_ccoeff(self.t_rspl, y_ccoeff, r)

            cur_coeff = np.log(self.maxrank/r)+1

            ax.plot((x_if-self.x_mean)*cur_coeff, (y_if-self.y_mean)*cur_coeff, color=mycolor, alpha=myalpha)

        ax.plot(self.x_val-self.x_mean, self.y_val-self.y_mean, color=self.dcolor)
    
        return(fig)



if __name__ == '__main__':

    data_in = 'input_fig/'
    dir_out = 'out/'

    # black background and one hue
    # assuming path closed, or at least closable
    fn = 'renard.csv' 
    my_sofourier = Sofourier(data_in + fn, delimiter=';', decimal=',')
    my_sofourier.set_color(main_color_hue=0.08)
    my_sofourier.set_maxrank(10)
    my_sofourier.set_sat(0, 1)
    my_sofourier.set_alpha(0.3, 1)
    myfig = my_sofourier.plot()
    myfig.savefig(dir_out + 'sofourier_renard.pdf')

    # white background and rainbow style

    fn_star = 'star.csv'
    my_sof = Sofourier(data_in + fn_star, delimiter=';', decimal=',')
    my_sof.set_color(color_mode='auto', data_color='black', background='white')
    my_sof.set_alpha(0.8, 1)
    my_sof.set_maxrank(10)
    myfig2 = my_sof.plot()
    myfig2.savefig(dir_out + 'sofourier_rainbowstar.pdf')