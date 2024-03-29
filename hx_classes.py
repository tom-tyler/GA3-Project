import numpy as np
from iapws._iapws import _Liquid
from ht.conv_tube_bank import baffle_leakage_Bell, unequal_baffle_spacing_Bell

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
                      plenum_length_1,
                      plenum_length_2,
                      baffle_gap,
                      baffle_type,
                      tube_layout,
                      shell_passes,
                      tube_bundle_diameter,
                      tube_passes,
                      baffle_spacing_in,
                      baffle_spacing_out,
                      design_year,
                      pump_year,
                      T_inh,
                      T_inc,
                      leakage = True,
                      name = None,
                      real_data = None,
                      ):

        #accuracy of calculations
        self.accuracy = 0.05
        self.T_increment = 5.0
        self.m_increment = 0.1

        #input heat exchanger parameters
        self.tube_number = tube_number #number of tubes
        self.baffle_number = baffle_number #number of baffles
        self.design_year = design_year
        self.pump_year = pump_year
        self.T_inh = T_inh
        self.T_inc = T_inc
        self.name = name
        self.real_data = real_data

        self.tube_layout = tube_layout #if tubes laid out in a triangular fashion, tube_layout = 't', if laid out in square fashion, tube_layout = 's'
        self.baffle_spacing_in = baffle_spacing_in
        self.baffle_spacing_out = baffle_spacing_out
        self.pitch = pitch #pitch between tubes
        self.tube_length = tube_length #length of copper tubes carrying fluid
        self.plenum_length_1 = plenum_length_1
        self.plenum_length_2 = plenum_length_2
        self.baffle_gap = baffle_gap #baffle gap
        self.leakage = leakage

        self.shell_passes = shell_passes #number of shell passes
        self.tube_passes = tube_passes


        #creating tube,shell,plate and baffle objects
        self.tube = pipe(6e-3,8e-3,0.20,3.5,self.tube_length)
        self.shell = pipe(64e-3,70e-3,0.650,0.5,self.tube_length)
        self.plenum1 = pipe(64e-3,70e-3,0.650,0.5,self.plenum_length_1)
        self.plenum2 = pipe(64e-3,70e-3,0.650,0.5,self.plenum_length_2)
        #self.supports = pipe(0, ) need to specify these
        self.plate = sheet(4.5e-3,6.375,self.shell.d_inner,circular = True)
        
        self.baffle_type = baffle_type
        if self.design_year in [2022,2020,2019]:
            nozzle_bore = 20e-3
        elif design_year in [2018,2017]:
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
            self.baffle_area = self.shell.c_area - self.segment_area - self.tube_number*(self.tube.c_area)

        elif baffle_type == 'across_b': #alternating circles and doughnuts, need to use a dictionary for baffle gap
            self.baffle_area_doughnut = self.shell.c_area - (self.shell.d_inner - 2*baffle_gap['doughnut'])**2*np.pi/4
            self.baffle_area_disk = (self.shell.d_inner - 2*baffle_gap['disk'])**2*np.pi/4
            self.baffle_area = np.mean(self.baffle_area_disk,self.baffle_area_doughnut) - self.tube_number*(self.tube.c_area)
        
        elif baffle_type == 'across_a':
            self.baffle_area = self.shell.c_area - self.tube_number*(self.tube.d_outer + baffle_gap)**2*np.pi/4

        else:
            print('specify correct baffle type please')
        self.baffle = sheet(1.5e-3,2.39,self.baffle_area)

        #derived quantities
        self.shell_length = self.tube_length - 2*self.plate.thickness
        self.nozzle_exit_area = (self.plenum_length_1)*(self.shell.d_inner) #'guess' at area which nozzle is exapnding into
        
        #if (self.baffle_spacing_in < 40e-3) and (self.baffle_spacing_out < 40e-3):
        if self.baffle_number == 0:
            print('please have at least one baffle lol')
        elif self.baffle_number == 1:
            print('please have at least one baffle lol')
            self.baffle_spacing_in, self.baffle_spacing_out, self.baffle_spacing = self.shell_length/2, self.shell_length/2, self.shell_length/2
        else:
            self.baffle_spacing = (self.shell_length - self.baffle_spacing_in - self.baffle_spacing_out)/(self.baffle_number - 1) #B, spacing between baffles
        #else:
            #self.baffle_spacing = self.shell_length/(self.baffle_number + 1)
        
        self.sigma = self.tube_number*self.tube.c_area/self.plate.area #sigma, used to find Ke and Kc
        self.sigma_nozzle = self.nozzle_c_area/(self.nozzle_exit_area) #sigma for nozzle
        self.convection_area = np.pi*self.tube.d_inner*(self.shell_length)*tube_number # total area of tube surface for convection
        self.baffle_cut = self.baffle_gap/self.shell.d_inner
        self.tube_bundle_diameter = tube_bundle_diameter
        
        #axisymmetric dividers for shell - multiple shell passes - only approximate to get mass
        if shell_passes > 1:
            self.divider_no = shell_passes
        else:
            self.divider_no = 0
        self.divider_area = self.shell_length*(self.shell.d_inner/2)
        self.divider = sheet(1.5e-3,2.39,self.divider_area)
        

        #for multi tube passes
        if tube_passes == 2:
            self.plen_dividers = 2
        elif tube_passes >= 3:
            self.plen_dividers = tube_passes+2
        else:
            self.plen_dividers = 0
        #approximate
        self.plen_divider_area = ((self.plenum_length_1 + self.plenum_length_2)/2)*(self.shell.d_inner/2)
        self.plen_divider = sheet(1.5e-3,2.39,self.plen_divider_area)

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
            self.pp = self.pitch*(np.sqrt(3)/2)
        elif self.tube_layout == 's':
            self.a1 = 0.321 
            self.a2 = -0.388
            self.a3 = 1.450
            self.a4 = 0.519
            self.theta = 0
            self.b1 = 0.0815
            self.b2 = 0.022
            self.b3 = 6.30
            self.b4 = 0.378
            self.area_adjustment_factor = 1
            self.pp = self.pitch
        elif self.tube_layout == 's_rot':
            self.a1 = 0.321 
            self.a2 = -0.388
            self.a3 = 1.450
            self.a4 = 0.519
            self.theta = 45
            self.b1 = 0.0815
            self.b2 = 0.022
            self.b3 = 6.30
            self.b4 = 0.378
            self.area_adjustment_factor = 1/(2**0.5)
            self.pp = self.pitch*(np.sqrt(2))

        #adjusted diameter parameters
        self.d_otl = tube_bundle_diameter
        self.D_ctl = self.d_otl - self.tube.d_outer

        #derived angle parameters
        arg_ctl = (self.shell.d_inner/self.D_ctl) * (1 - 2*self.baffle_cut)
        if abs(arg_ctl) > 1:
            arg_ctl = 0.99
        self.theta_ctl = 2*np.arccos(arg_ctl)/self.shell_passes
        self.theta_ds = 2*np.arccos(1 - 2*self.baffle_cut)/self.shell_passes

        #F factors
        self.Fw = (1/(2*np.pi))*(self.theta_ctl - np.sin(self.theta_ctl))
        self.Fc = 1 - 2*self.Fw

        #clearances
        self.delta_tb = 0
        self.delta_sb = 2.5e-4

        #effective areas
        self.Sm = self.baffle_spacing * ((self.shell.d_inner - self.d_otl) + ((self.d_otl - self.tube.d_outer)*(self.pitch - self.tube.d_outer))/(self.pitch*self.area_adjustment_factor))/self.shell_passes
        self.Sb = self.baffle_spacing * (self.shell.d_inner - self.d_otl - self.tube.d_outer/2)/self.shell_passes
        self.Ssb = self.shell.d_inner * self.delta_sb* (np.pi - 0.5*self.theta_ds)
        self.Stb = (np.pi/4) * ((self.tube.d_outer + self.delta_tb)**2 - self.tube.d_outer**2) * self.tube_number * (1 - self.Fw)
        self.A_shell = self.Sm

        #derived area ratios
        self.rs = self.Ssb / (self.Ssb + self.Stb)
        self.rl = (self.Ssb + self.Stb) / self.Sm
        self.p = 0.8 - 0.15*(1 + self.rs)
    
        #pressure drop variables
        #region
        self.Nc = round(self.shell.d_inner * (1 - 2 * self.baffle_cut) / (self.pp))
        self.Ncw = round(0.8*(baffle_gap) / (self.pp))
        if self.Ncw < 0:
            self.Ncw = 0
        if self.Nc < 1:
            self.Nc = 1
        if leakage == True:
            self.Rl = np.exp(-1.33*(1+self.rs)*(self.rl**self.p))
        else:
            self.Rl = 1 #assume no leakage
        self.Rb = np.exp(-3.7*(self.Sb/self.Sm))
        if shell_passes == 2:
            self.Rs = ((self.baffle_spacing/self.baffle_spacing_in)**(1.8))
        elif shell_passes == 1:
            self.Rs = 0.5*((self.baffle_spacing/self.baffle_spacing_in)**(1.8) + (self.baffle_spacing/self.baffle_spacing_out)**(1.8))
        #endregion

        #thermal design variables
        #region
        self.Jc = 0.55 + 0.72*self.Fc
        if leakage == True:
            self.Jl = baffle_leakage_Bell(self.Ssb,self.Stb,self.Sm)
        else:
            self.Jl = 1 #assume no leakage
        self.Jb = np.exp(-1.25*(self.Sb/self.Sm))
        self.Jr = 1 #=1 due to high Re
        self.Js = unequal_baffle_spacing_Bell(self.baffle_number,self.baffle_spacing,self.baffle_spacing_in,self.baffle_spacing_out)
        #endregion       
        

    def total_mass(self):
        #calculate total mass of heat exchanger
        self.total_nozzle_mass = self.no_nozzles*self.nozzle_mass
        self.total_tube_mass = self.tube_number*self.tube.mass
        self.total_shell_mass = self.shell.mass + self.plenum1.mass + self.plenum2.mass
        self.total_baffle_mass = self.baffle_number*self.baffle.mass
        self.total_plate_mass = self.plate_number*self.plate.mass
        self.total_divider_mass = self.divider_no*self.divider.mass + self.plen_dividers*self.plen_divider.mass

        return self.total_baffle_mass + self.total_tube_mass + self.total_nozzle_mass + self.total_plate_mass + self.total_shell_mass + self.total_divider_mass
