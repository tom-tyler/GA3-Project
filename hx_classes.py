import numpy as np
from iapws._iapws import _Liquid
from ht.conv_tube_bank import baffle_correction_Bell, baffle_leakage_Bell , bundle_bypassing_Bell, laminar_correction_Bell, unequal_baffle_spacing_Bell
from ht.conv_tube_bank import dP_Kern, dP_Zukauskas

class water:
    def __init__(self,Tin,Tout,Pin = 1e5, Pout = 1e5):
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
                      crossflow_tube_fraction,
                      bypass_area_fraction,
                      seal_strips,
                      crossflow_rows,
                      tube_bundle_diameter,
                      tube_passes,
                      year,
                      T_inh,
                      T_inc,
                      co_counter='counter',
                      approximate_glue_mass=0
                      ):

        #input heat exchanger parameters
        self.tube_number = tube_number #number of tubes
        self.baffle_number = baffle_number #number of baffles
        self.year = year
        self.T_inh = T_inh
        self.T_inc = T_inc

        self.tube_layout = tube_layout #if tubes laid out in a triangular fashion, tube_layout = 't', if laid out in square fashion, tube_layout = 's'
        self.crossflow_tube_fraction = crossflow_tube_fraction #fraction of tubes which are between baffle tips and not in the window
        self.bypass_area_fraction = bypass_area_fraction #Fraction of the crossflow area which is not blocked by a baffle or anything else and available for bypassing,
        self.seal_strips = seal_strips #Number of seal strips per side of a baffle added to prevent bypassing
        self.crossflow_rows = crossflow_rows #The number of tube rows in the crossflow of the baffle

        self.pitch = pitch #pitch between tubes
        self.tube_length = tube_length #length of copper tubes carrying fluid
        self.shell_length = shell_length #length of shell
        self.baffle_gap = baffle_gap #baffle gap

        self.shell_passes = shell_passes #number of shell passes
        self.tube_passes = tube_passes

        self.co_counter = co_counter #counter or co flow. for counter = 'counter', for co = 'co'
        self.approximate_glue_mass = approximate_glue_mass #approximate mass of glue (may move this to fixed parameters)


        #creating tube,shell,plate and baffle objects
        self.tube = pipe(6e-3,8e-3,0.20,3.5,self.tube_length)
        self.shell = pipe(64e-3,70e-3,0.650,0.5,self.shell_length)
        #self.supports = pipe(0, ) need to specify these
        self.plate = sheet(4.5e-3,6.375,self.shell.d_inner,circular = True)
        
        self.baffle_type = baffle_type
        if self.year in [2022,2020,2019]:
            nozzle_bore = 20e-3
        elif year in [2018,2017]:
            nozzle_bore = 19e-3
        self.nozzle_bore = nozzle_bore #bore size of nozzle

        #fixed heat exchanger parameters
        self.nozzle_c_area = (np.pi*self.nozzle_bore**2)/4
        self.no_nozzles = 4 #number of nozzles
        self.nozzle_mass = 0.025 #mass of nozzles
        self.plate_number = 4 #number of plates

        if baffle_type == 'across_c': #normal case
            r = self.shell.d_inner/2
            self.segment_area = r**2 * np.arccos((r - baffle_gap)/r) - (r - baffle_gap)*(2*r*baffle_gap - baffle_gap**2)**0.5
            self.baffle_area = self.shell.c_area - self.segment_area

        elif baffle_type == 'across_b': #alternating circles and doughnuts, need to use a dictionary for baffle gap
            self.baffle_area_doughnut = self.shell.c_area - (self.shell.d_inner - 2*baffle_gap['doughnut'])**2*np.pi/4
            self.baffle_area_disk = (self.shell.d_inner - 2*baffle_gap['disk'])**2*np.pi/4
            self.baffle_area = np.mean(self.baffle_area_disk,self.baffle_area_doughnut)
        
        elif baffle_type == 'across_a':
            self.baffle_area = self.shell.c_area - self.tube_number*(self.tube.d_outer + baffle_gap)**2*np.pi/4

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
        self.tube_length_in_shell = self.tube_length - 2*self.plate.thickness
        self.convection_area = np.pi*self.tube.d_inner*(self.tube_length_in_shell)*tube_number # total area of tube surface for convection
        self.baffle_cut = self.baffle_gap/self.shell.d_inner
        self.tube_bundle_diameter = tube_bundle_diameter
        self.d_otl = tube_bundle_diameter
        self.D_ctl = self.d_otl - self.tube.d_outer
        
        #axisymmetric dividers
        if shell_passes > 1:
            self.divider_no = shell_passes
        else:
            self.divider_no = 0
        self.divider_area = self.tube_length_in_shell*(self.shell.d_inner/2)
        self.divider = sheet(1.5e-3,2.39,self.divider_area)

        #bell delaware coefficients
        if self.tube_layout == 't':
            self.a1 = 0.321 
            self.a2 = -0.388
            self.a3 = 1.450
            self.a4 = 0.519
            self.theta = 30
            self.b1 = 0.372
            self.b2 = -0.123
            self.b3 = 7
            self.b4 = 0.5
            self.area_adjustment_factor = 1
        elif self.tube_layout == 's':
            self.theta = 0
            self.b1 = 0.0815
            self.b2 = 0.022
            self.b3 = 6.30
            self.b4 = 0.378
            self.area_adjustment_factor = 1
        elif self.tube_layout == 's_rot':
            self.theta = 45
            self.area_adjustment_factor = 1/(2**0.5)

        self.theta_ds = 2*np.arccos(1 - 2*self.baffle_cut)
        self.theta_ctl = 2*np.arccos((self.shell.d_inner/self.D_ctl) * (1 - 2*self.baffle_cut))
        self.Fw = (1/(2*np.pi))*(self.theta_ctl - np.sin(self.theta_ctl))
        self.Fc = 1 - 2*self.Fw
        self.A_tubes = self.tube_number * self.Fw * (np.pi*self.tube.d_outer**2/4)

        #effective areas
        self.Sm = self.baffle_spacing * ((self.shell.d_inner - self.d_otl) + ((self.d_otl - self.tube.d_outer)*(self.pitch - self.tube.d_outer))/(self.pitch*self.area_adjustment_factor))
        self.Swg = (1/8) * self.shell.d_inner**2 * (self.theta_ds - np.sin(self.theta_ds)) 
        self.Sw = self.Swg - self.A_tubes
        self.Sb = self.baffle_spacing * (self.shell.d_inner - self.d_otl)
        self.Nc = self.shell.d_inner * (1 - 2 * self.baffle_cut) / (self.pitch * np.cos(self.theta))
        self.Ncw = (0.8 * self.baffle_cut * self.shell.d_inner) / (self.pitch * np.cos(self.theta))

        #clearances
        self.delta_tb = 0
        self.delta_sb = 0.8e-3 + 0.002 * self.shell.d_inner

        #self.Ssb = self.shell.d_inner * self.delta_sb* (np.pi - 0.5*self.theta_ds)
        #self.Stb = np.pi * self.tube.d_outer * self.delta_tb * self.tube_number * (1+self.Fc)
        self.Ssb = 0
        self.Stb = 0

        #self.rs = self.Ssb / (self.Ssb + self.Stb)
        self.rs = 1
        #self.rl = (self.Ssb + self.Stb) / self.Sm
        self.rl = 0
        self.Rs = 1 #uniform baffle spacing
        self.Rb = np.exp(-3.7*(self.Sb/self.Sm))
        self.p = 0.8 - 0.15*(1 + self.rs)
        #self.Rl = np.exp(-1.33*(1+self.rs)*(self.rl**self.p))
        self.Rl = 1 # no leakage

        #self.Jc = baffle_correction_Bell(self.crossflow_tube_fraction, method = 'chebyshev')
        self.Jc = 0.55 + 0.72*self.Fc
        #self.Jl = baffle_leakage_Bell(self.Ssb,self.Stb,self.Sm)
        self.Jl = 1 #no leakage
        #self.Jb = bundle_bypassing_Bell(self.bypass_area_fraction,self.seal_strips,self.crossflow_rows)
        self.Jb = 1
        self.Jr = 1 #=1 due to high Re
        self.Js = 1 #=1 due to even baffle spacing

    def total_mass(self):
        #calculate total mass of heat exchanger
        self.total_nozzle_mass = self.no_nozzles*self.nozzle_mass
        self.total_tube_mass = self.tube_number*self.tube.mass
        self.total_shell_mass = self.shell.mass
        self.total_baffle_mass = self.baffle_number*self.baffle.mass
        self.total_plate_mass = self.plate_number*self.plate.mass
        self.total_divider_mass = self.divider_no*self.divider.mass
        #self.total_mass_supports = 

        return self.total_baffle_mass + self.total_tube_mass + self.total_nozzle_mass + self.total_plate_mass + self.total_shell_mass + self.total_divider_mass
