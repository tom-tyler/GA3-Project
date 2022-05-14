import numpy as np
from iapws._iapws import _Liquid
import hx_functions as hxf

class water:
    def __init__(self,Tin,Tout,Pin,Pout):
        self.Tin = Tin
        self.Tout = Tout
        self.Pin = Pin
        self.Pout = Pout
        self.avg_T = (Tin + Tout)/2
        self.avg_P = (Pin + Pout)/2

        w = _Liquid(self.avg_T + 273.15, self.avg_P / 1e6) #needs to be in degrees C and Pa
        self.cp = w['cp']
        self.k = w['k']
        self.mu = w['mu']
        self.Pr = w['cp']*w['mu']/w['k']
        self.rho = w['rho']

class pipe:
    def __init__(self,d_inner,d_outer,rho_l,max_l,l):
        self.d_inner = d_inner
        self.d_outer = d_outer
        self.rho_l = rho_l
        self.max_l = max_l
        self.l = l
        self.mass = self.l*self.rho_l
        self.c_area = (np.pi*d_inner**2)/4

class sheet:
    def __init__(self,thickness,rho_A,geometry,circular = False):
        self.thickness = thickness
        self.rho_A = rho_A

        if circular == True:
            self.diameter = geometry
            self.area = (np.pi*self.diameter**2)/4
        else:
            self.area = geometry
        self.mass = self.area*self.rho_A

class HX:
    def __init__(self,tube_number,baffle_number,pitch,tube_length,shell_length,baffle_area,tube_layout,shell_passes,approximate_glue_mass=0):

        #input heat exchanger parameters
        self.tube_number = tube_number #number of tubes
        self.baffle_number = baffle_number #number of baffles
        self.pitch = pitch #pitch between tubes
        self.tube_layout = tube_layout #if tubes laid out in a triangular fashion, tube_layout = 't', if laid out in square fashion, tube_layout = 's'
        self.tube_length = tube_length #length of copper tubes carrying fluid
        self.shell_length = shell_length #length of shell
        self.baffle_area = baffle_area #cross-sectional area of baffle
        self.shell_passes = shell_passes
        self.approximate_glue_mass = approximate_glue_mass #approximate mass of glue (may move this to fixed parameters)

        #fixed heat exchanger parameters
        self.nozzle_bore = 20e-3 #bore size of nozzle
        self.nozzle_c_area = (np.pi*self.nozzle_bore**2)/4
        self.no_nozzles = 4 #number of nozzles
        self.nozzle_mass = 0.025 #mass of nozzles
        self.plate_number = 4 #number of plates

        #creating tube,shell,plate and baffle objects
        self.tube = pipe(6e-3,8e-3,0.20,3.5,self.tube_length)
        self.shell = pipe(64e-3,70e-3,0.650,0.5,self.shell_length)
        self.plate = sheet(4.5e-3,6.375,self.shell.d_inner,circular = True)
        self.baffle = sheet(1.5e-3,2.39,self.baffle_area)

        #derived quantities
        self.baffle_spacing = self.tube_length/(1+self.baffle_number) #B, spacing between baffles
        self.sigma = self.tube_number*self.tube.c_area/self.plate.area #sigma, used to find Ke and Kc
        self.A_shell = self.shell.d_inner*(self.pitch-self.tube.d_outer)*self.baffle_spacing/self.pitch #area through which shell fluid can flow

        
    def total_mass(self):
        #calculate total mass of heat exchanger
        self.total_nozzle_mass = self.no_nozzles*self.nozzle_mass
        self.total_tube_mass = self.tube_number*self.tube.mass
        self.total_shell_mass = self.shell.mass
        self.total_baffle_mass = self.baffle_number*self.baffle.mass
        self.total_plate_mass = self.plate_number*self.plate.mass

        return self.total_baffle_mass + self.total_tube_mass + self.total_nozzle_mass + self.total_plate_mass + self.total_shell_mass

