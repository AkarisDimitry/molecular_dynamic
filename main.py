#-*- coding: utf-8 -*-
from particles.Particles import *
from interaction.Interaction import * 
from integrator.Integrator import *
from thermostat.Thermostat import *
from box.Box import *
#from interaction.descriptor.Descriptors import *

import time
import matplotlib
import matplotlib.pyplot as plt
from  mpl_toolkits.mplot3d import axes3d

import numpy as np

def simulation(particles, interaction_parameters, integration, boxdata, thermostat, steps, save=True):
   
   # ---- PARTICLE ---- #
   # make particle system
   sistema = Base() 
   sistema.load(particles)
   sistema.alocate_memory(steps=steps)
   
   # ---- INTERACTION ---- #
   # make interaction object 
   interaction = interaction_config( interaction_parameters )

   # ---- INTEGRATION ---- #
   # make integration object
   if integration['type'] == 'Velverlet':
      integrator = VelVerlet( integration['dt'] )

   # ---- BOX ---- #
   # make periodic conditions
   box = Box( x0=boxdata['min'], xf=boxdata['max'], t=boxdata['type'] )

   # ---- THERMOSTAT ---- #
   # set the thermostat
   thermostat = thermostat 

   # ---- STEPS ---- #
   # set total steps
   step, steps = np.zeros(1, dtype=np.int), steps

   # ---- SAVE ---- #
   # save
   save = True

   while step < steps:
      # --- measure time --- #
      t1 = time.clock()

      # --- Integrator first step --- #    
      t2 = time.clock()
      sistema.x, sistema.v = integrator.first_step(sistema.x, sistema.v, sistema.f, sistema.constrains)
      #print('Primer integrador', time.clock()-t2)

      # --- Forces --- #
      t3 = time.clock()
      sistema.f, sistema.E = interaction.Cforces(sistema.x, sistema.v, sistema.charge, sistema.t, box.xf)
      #print('Forces', time.clock()-t3)

      # --- Store data --- #
      t4 = time.clock()
      if save:  sistema.store_step(step=step)
      #print('Sistema', time.clock()-t4)

      # --- Integrator last step --- #
      t5 = time.clock()
      sistema.x, sistema.v = integrator.last_step(sistema.x, sistema.v, sistema.f, sistema.constrains)    
      #print('Integrador', time.clock()-t5)
      
      # --- BOX aprox --- #
      t6 = time.clock()
      sistema.x, sistema.v = box.C_wrap_boundary(sistema.x, sistema.v)
      #print('Box', time.clock()-t6)

      sistema.v = thermostat.friction(v=sistema.v)

      step += 1
      if step%100 == 0: print(step)
      #print(step, ' *** Total time *** ', t1 - time.clock())

   sistema.save_OUTFILE( particles['OUTFILE'] )

   if True:
      plt.plot(sistema.store_data['Ec'] )
      plt.plot(sistema.store_data['Ep'] )
      plt.plot(np.array(sistema.store_data['Ec'])+np.array(sistema.store_data['Ep']) )
      plt.show()



# --- INCAR --- #
particles = {              'init':              'random', 
                           'n':                 15, 
                           'filename':          'POSFILE',
                           'dimentionality':    2,
                           'OUTFILE':           'CONFILE',
                           'position_radius':   100}

particles = {              'init':              'material', 
                           'RD':                 {'a1':[1.3,0,0], 'a2':[0,1.3,0], 'a3':[0,0,1]},
                           'rep':                {'x':5, 'y':5, 'z':5},
                           'OUTFILE':           'CONFILE',}

interaction_parameters =   {  
                              'type':  'flocking', 
                              'rcut':  10.5, 
                              'coef01':  [0,1,0],
                              'coef02':  [0,-100,0],
                              'coef03':  [0,-3,0], 
                              'shift_style':'Displace' }


interaction_parameters =   {  'type':  'Coulomb', 
                              'rcut':  5.5, 
                              'eps':   1.0, 
                              'sigma': 1.0, 
                              'shift_style':'Displace' }

interaction_parameters =   {  'type':  'LJ', 
                              'rcut':  5.5, 
                              'eps':   -2.0, 
                              'sigma': 1.0, 
                              'shift_style':'Displace' }

integration= {
               'type':  'Velverlet', 
               'dt':    0.05} 

boxdata = { 
            'type': 'Periodic',
            'min': [-0,-0,-0], 
            'max': [ 1.3*6, 1.3*6, 300], }

thermostat = FriccionCoef(coef=0.9, dt=integration['dt'])
steps = 2000

simulation( particles=particles, 
            interaction_parameters=interaction_parameters, 
            integration=integration, 
            boxdata=boxdata, 
            thermostat=thermostat, 
            steps=steps)
'''
# todo 
(1) clean interacion code
(2) add INFILE
(3) share distance matrix in memory pointer

'''