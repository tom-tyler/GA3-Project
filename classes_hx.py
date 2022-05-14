class HX:
    def __init__(self):
        self.tube_number = [] 
        self.baffle_number = []
        self.pitch = []
        self.nozzle_diameter = 20e-3
        self.tube_outer_diameter = 8e-3
        self.tube_inner_diameter = 6e-3
        self.tube_weight_per_length = 0.20
        
        
    # This is a method
    def full_name(self):
        "Return full name (forename surname)"
        return self.forename + " " + self.surname

    # This is a method
    def surname_forename(self, sep=","):
        "Return 'surname, forename', with option to specify separator"
        return self.surname + sep + " " + self.forename