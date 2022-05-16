"""
Display module
"""
import numpy as np
import pygame
from sys import exit
from pygame.locals import *
class Display(object):
    """
    Base display class
    """
    def __init__(self, metadata=[0,'Unknow','Unknow','None'], resolution=(600,600), FPS = 60, flags=pygame.DOUBLEBUF):
        self.resolution = resolution
        self.flags = flags 
        self.screen = pygame.display.set_mode( (self.resolution[0], self.resolution[1]), self.flags )
        self.FPS = FPS 
        self.zoom = [-260,-300,10]
        self.time = 0
        self.atomic_data = []
        self.graph_data = []

        self.meta_data = [
                'Particle number: '+str(metadata[0]),
                'Inicialization: '+str(metadata[1]),
                'potencial: '+str(metadata[2]),
                'Other: '+str(metadata[3])] 
            # (0) numero de particulas 
            # (1) inicializacion
            # (2) potencial
            # (3) other info 

    def update(self):
        data = self.atomic_data
        timer = pygame.time.Clock()
        while True:
            if int(self.time)+1 < data.shape[2]: self.time += 0.5
            else: self.time = 0
            self.screen.fill((0,0,0))
            t = int(self.time)
            for m in range(4):
              Xdata = data[:,:,:,m]
              for n in range(data.shape[1]):
                x = Xdata[0][n][t]*50 / Xdata[2][n][t] * self.zoom[2] + self.zoom[0] 
                y = Xdata[1][n][t]*50 / Xdata[2][n][t] * self.zoom[2] + self.zoom[1]
                r = 5
                
                color = (30+220*float(Xdata[2,n,t]-np.min(Xdata[2,:,:]))/np.max(Xdata[2,:,:]-np.min(Xdata[2,:,:])), 0, m*50)
                 
                pygame.draw.circle(self.screen, color, ( int(x), int(y)), int(r) )
            
            if t > 2:
                for n in self.graph_data:
                   f1, f2 = float(self.resolution[0]-200)/n.shape[0], 200.0/(np.max(n)-np.min(n) )
                   m = (n[:t]-np.min(n))*f2
                   pygame.draw.lines(self.screen, (255,255,255), 0, zip(
                       np.arange(t)*f1+100,
                       -m+self.resolution[1]-100), 2 )

            
            
            
            for i, n in enumerate(self.meta_data):
                self.screen.blit(pygame.font.Font(None, 20 ).render(n, 0, (0,255,0)),[50,400+i*20])
            self.screen.blit(pygame.font.Font(None, 20 ).render('Time: ' + str(self.time), 0, (0,255,0)),[400,400])           
            self.screen.blit(pygame.font.Font(None, 20 ).render('Step: ' + str(int(self.time)), 0, (0,255,0)),[400,400+20])
            self.screen.blit(pygame.font.Font(None, 20 ).render('Zoom data: ' + str(int(self.zoom[0]))+', '+str(int(self.zoom[1])) +', '+str(int(self.zoom[2])), 0, (0,255,0)),[400,400+40])
            
            
                
            timer.tick(self.FPS)
            pygame.display.update()
    
    def transform(x,y,z,v1,v2):
        a = 1

    def change_zoom(self, cursor, f ):
        if cursor.key[0] == 1:
            if self.zoom[2] < 100000:
                self.zoom[2] *= f
                self.zoom[1] = -float((( cursor.top - self.zoom[1]) * f) - cursor.top)
                self.zoom[0] = -float((( cursor.left - self.zoom[0]) * f) - cursor.left)
        if cursor.key[0] == -1:
            self.zoom[2] = float(self.zoom[2])/f 
            if self.zoom[2] <= 0: self.zoom[2] = 1
            self.zoom[1] = -float(((cursor.top - self.zoom[1]) / f) - cursor.top)
            self.zoom[0] = -float(((cursor.left - self.zoom[0]) / f) - cursor.left) 
        
    def keys_press(self, cursor ):
        for event in pygame.event.get():
            if event.type == pygame.MOUSEBUTTONDOWN:
                cursor.update()
                cursor.mouse[event.button] = 1
                cursor.mouse_clickposition = [[cursor.left, self.zoom[0]],[cursor.top, self.zoom[1]]]
            if event.type == pygame.MOUSEMOTION:
                cursor.update()
                if cursor.mouse[1] == 1:
                    self.zoom[0] = cursor.mouse_clickposition[0][1] - (cursor.mouse_clickposition[0][0]- cursor.left)
                    self.zoom[1] = cursor.mouse_clickposition[1][1] - (cursor.mouse_clickposition[1][0]- cursor.top)

            if event.type == pygame.MOUSEBUTTONUP:
                cursor.mouse[event.button] = 0 

            if event.type == QUIT:
                pygame.quit()
               
            if event.type == pygame.KEYDOWN:
                 if event.key == pygame.K_q:
                     pygame.quit()
            
                 if event.key == pygame.K_DOWN:
                     cursor.key[0] = -1

                 if event.key == pygame.K_UP:
                     cursor.key[0] = 1

            if event.type == pygame.KEYUP:
                 if event.key == pygame.K_DOWN:
                     cursor.key[0] = 0

                 if event.key == pygame.K_UP:
                     cursor.key[0] = 0




class Cursor(object):
    def __init__(self):
        self.left, self.top = pygame.mouse.get_pos()
        self.key = [0,0,0,0,0]
        self.mouse = [0,0,0,0,0,0]
        self.mouse_clickposition = []
    def update(self):
        self.left, self.top = pygame.mouse.get_pos()




















