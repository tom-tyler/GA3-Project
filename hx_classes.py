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
        self.cp = w['cp'] *1000
        self.k = w['k']
        self.mu = w['mu']
        self.Pr = self.cp*self.mu/self.k
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
    def __init__(self,tube_number,
                      baffle_number,
                      pitch,
                      tube_length,
                      shell_length,
                      baffle_gap,
                      baffle_type,
                      tube_layout,
                      shell_passes,
                      tube_passes,
                      nozzle_bore,
                      co_counter='counter',
                      approximate_glue_mass=0
                      ):

        #input heat exchanger parameters
        self.tube_number = tube_number #number of tubes
        self.baffle_number = baffle_number #number of baffles

        self.tube_layout = tube_layout #if tubes laid out in a triangular fashion, tube_layout = 't', if laid out in square fashion, tube_layout = 's'
        
        self.pitch = pitch #pitch between tubes
        self.tube_length = tube_length #length of copper tubes carrying fluid
        self.shell_length = shell_length #length of shell
        self.baffle_area = baffle_gap #cross-sectional area of baffle

        self.F = 1 #this is to be changed later

        self.tube_passes = tube_passes #number of tube passes
        self.shell_passes = shell_passes #number of shell passes

        self.co_counter = co_counter #counter or co flow. for counter = 'counter', for co = 'co'
        self.approximate_glue_mass = approximate_glue_mass #approximate mass of glue (may move this to fixed parameters)
        self.nozzle_bore = nozzle_bore #bore size of nozzle

        #fixed heat exchanger parameters
        self.nozzle_c_area = (np.pi*self.nozzle_bore**2)/4
        self.no_nozzles = 4 #number of nozzles
        self.nozzle_mass = 0.025 #mass of nozzles
        self.plate_number = 4 #number of plates

        #creating tube,shell,plate and baffle objects
        self.tube = pipe(6e-3,8e-3,0.20,3.5,self.tube_length)
        self.shell = pipe(64e-3,70e-3,0.650,0.5,self.shell_length)
        #self.supports = pipe(0, ) need to specify these
        self.plate = sheet(4.5e-3,6.375,self.shell.d_inner,circular = True)
        
        self.baffle_type = baffle_type #across shell denoted by 'across'. Along shell denoted by 'along', along shell, but axisymmetric denoted by 'axial'
        
        if baffle_type == 'across_c': #normal case
            r = self.shell.d_inner/2
            self.segment_area = r**2 * np.arccos((r - baffle_gap)/r) - (r - baffle_gap)*(2*r*baffle_gap - baffle_gap**2)**0.5
            self.baffle_area = self.shell.c_area - self.segment_area

        elif baffle_type == 'across_b': #alternating circles and doughnuts, need to use a dictionary for baffle gap
            self.baffle_area_doughnut = self.shell.c_area - (self.shell.d_inner - 2*baffle_gap['doughnut'])**2*np.pi/4
            self.baffle_area_disk = (self.shell.d_inner - 2*baffle_gap['disk'])**2*np.pi/4
            self.baffle_area = np.mean(self.baffle_area_disk,self.baffle_area_doughnut)
        
        elif baffle_type == 'across_a':
            self.baffel_area = self.shell.c_area - self.tube_number*(self.tube.d_outer + baffle_gap)**2*np.pi/4

        elif baffle_type == 'along':
            self.baffle_area = 0 #need to do this
        elif baffle_type == 'axial':
            self.baffle_area = 0 #need to do this
        else:
            print('specify correct baffle type please')
        self.baffle = sheet(1.5e-3,2.39,self.baffle_area)


        #derived quantities
        self.shell_end_length = (shell_length-tube_length)/2 #length of plenum
        self.nozzle_exit_area = (self.shell_end_length)*(self.nozzle_bore + self.plate.diameter)/2 #'guess' at area which nozzle is exapnding into
        self.baffle_spacing = self.tube_length/(1 + self.baffle_number) #B, spacing between baffles
        self.sigma = self.tube_number*self.tube.c_area/self.plate.area #sigma, used to find Ke and Kc
        self.sigma_nozzle = self.nozzle_c_area/(self.nozzle_exit_area) #sigma for nozzle
        self.A_shell = self.shell.d_inner*(self.pitch - self.tube.d_outer)*self.baffle_spacing/self.pitch #area through which shell fluid can flow
        self.convection_area = np.pi*self.tube.d_inner*(self.tube_length - 2*self.plate.thickness)*tube_number # total area of tube surface for convection
        
    def total_mass(self):
        #calculate total mass of heat exchanger
        self.total_nozzle_mass = self.no_nozzles*self.nozzle_mass
        self.total_tube_mass = self.tube_number*self.tube.mass
        self.total_shell_mass = self.shell.mass
        self.total_baffle_mass = self.baffle_number*self.baffle.mass
        self.total_plate_mass = self.plate_number*self.plate.mass
        #self.total_mass_supports = 

        return self.total_baffle_mass + self.total_tube_mass + self.total_nozzle_mass + self.total_plate_mass + self.total_shell_mass

