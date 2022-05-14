import numpy as np

class pipe:
    def __init__(self,d_inner,d_outer,rho_l,max_l,l):
        self.d_inner = d_inner
        self.d_outer = d_outer
        self.rho_l = rho_l
        self.max_l = max_l
        self.l = l
        self.mass = self.l*self.rho_l
        self.c_area = np.pi*self.d_inner**2/4

class sheet:
    def __init__(self,thickness,rho_A,A):
        self.thickness = thickness
        self.rho_A = rho_A
        self.A = A
        self.mass = self.A*self.rho_A

class HX:
    def __init__(self):

        #input heat exchanger parameters
        self.tube_number = [] 
        self.baffle_number = []
        self.plate_number = []
        self.pitch = []
        self.tube_length = []
        self.shell_length = []
        self.baffle_area = []
        self.plate_area = []
        self.approximate_glue_mass = []

        #fixed heat exchanger parameters
        self.tube = pipe(6e-3,8e-3,0.20,3.5,self.tube_length)
        self.shell = pipe(64e-3,70e-3,0.650,0.5,self.shell_length)
        self.nozzle_bore = 20e-3
        self.no_nozzles = 4
        self.nozzle_mass = 0.025
        self.plate = sheet(4.5e-3,6.375,self.baffle_area)
        self.baffle = sheet(1.5e-3,2.39,self.plate_area)

        
    def total_mass(self):
        #calculate total mass of heat exchanger
        self.total_nozzle_mass = self.no_nozzles*self.nozzle_mass
        self.total_tube_mass = self.tube_number*self.tube.mass
        self.total_shell_mass = self.shell.mass
        self.total_baffle_mass = self.baffle_number*self.baffle.mass
        self.total_plate_mass = self.plate_number*self.plate.mass

        return self.total_baffle_mass + self.total_tube_mass + self.total_nozzle_mass + self.total_plate_mass + self.total_shell_mass
    
    def sigma(self):

        return self.tube_number*self.tube.c_area/self.plate.A


    # This is a method
    def surname_forename(self, sep=","):
        "Return 'surname, forename', with option to specify separator"
        return self.surname + sep + " " + self.forename