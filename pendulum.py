#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 17:53:29 2020

@author: Isaac Ayala

See LICENSE for details.

Sample code to simulate a simple pendulum in Python


--------------
Explanation
--------------

1) Import libraries needed for the simulation.
2) Define custom functions.
3) Execute main script.


"""

from matplotlib import pyplot as plt
from  matplotlib import animation as animation
import numpy as np
from scipy import integrate

from celluloid import Camera



def system_dynamics(t, x, params,):
    """
    Parameters
    ----------
    x0 : State vector
    t : Current time step
    params : Simulation parameters

    Returns
    -------
    dx : State vector dynamics for time step integration

    """
    
    
    # Extract state variables and parameters
    # Python starts counting with 0 like any sane programming language!
    x1, x2, = x
    
    
    # Params is a dictionary with key and value pairs
    m = params['mass']
    g = params['gravity']
    l = params['length']
    k = params['friction']
    
    # Solve system dynamics for the current time step
    dot_x1 = x2
    dot_x2 = - (g/l) * np.sin(x1) - (k/m) * x2
    
    # Store each state variable into array
    dx = np.array([dot_x1, dot_x2])
    
    return dx

def plot_solution(sol, params):
    # Set figure size in inches
    fig = plt.figure(figsize=(12, 6))
    
    ax = fig.add_subplot(1, 2, 1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

    linear_plot = fig.add_subplot(2, 2, 2)
    linear_plot.set_xlabel(r'$t$')
    linear_plot.set_ylabel(r'$\chi$')

    angular_plot = fig.add_subplot(2, 2, 4)
    angular_plot.set_xlabel(r'$t$')
    angular_plot.set_ylabel(r'$\xi$')
    
    
    l = params['length']   
    dt = params['step']
    x = l * np.sin(sol.y[0])
    y = -l * np.cos(sol.y[0])
    
    ax.plot(x, y, '-', label = r'Linear')
    ax.set_xlim([-1.1*l, 1.1*l])
    ax.set_ylim([-1.1*l, 1.1*l])
    ax.legend()
    
    linear_plot.plot(sol.t, x, '-', label = r'$x$')
    linear_plot.plot(sol.t, y, '-', label = r'$y$')
    linear_plot.legend()
    
    
    
    angular_plot.plot(sol.t, sol.y[0], '-', label = r'$\theta$')
    angular_plot.plot(sol.t, sol.y[1], '-', label = r'$\dot \theta$')
    angular_plot.legend()
    
    fig.tight_layout()
    
    # Animation
    
    fig2 = plt.figure()
    anim = fig2.add_subplot(111, autoscale_on=False, xlim=(-1.1*l, 1.1*l), ylim=(-1.1*l, 1.1*l))
    anim.grid()

    line, = anim.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = anim.text(0.05, 0.9, '', transform=anim.transAxes)
    
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text


    def animate(i):
        thisx = [0, x[i]]
        thisy = [0, y[i]]
    
        line.set_data(thisx, thisy)
        time_text.set_text(time_template % (i*dt))
        return line, time_text
    
    ani = animation.FuncAnimation(fig2, animate, np.arange(1, len(sol.t)),
                              interval=25, blit=True, init_func=init)
    
    plt.show()
    

if __name__ == '__main__':
    # Close figures
    plt.close("all")
    
    # Simulation parameters
    dt = 0.01
    max_time = 10 
    t_span = [0, max_time] 
    t = np.arange(0, max_time, dt) 
    
    params = dict(
                    gravity = 9.81, # m/s^2
                    mass = 0.3, # kg
                    length = 0.15, # m
                    friction = 0.35,
                    step = dt,
                    )
    
    
    # Initial vector state
    x0 = np.array([np.pi/2, 0])
    
    sol = integrate.solve_ivp(
                                system_dynamics,
                                t_span,
                                x0,
                                method='RK45',
                                t_eval = t,
                                args = (params, )
                                )
    plot_solution(sol, params)
    
    
