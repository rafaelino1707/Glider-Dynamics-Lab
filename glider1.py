import math

class Glider:
    def __init__(self, mass, chord, wingspan, oswald_efficiency, cd0):
        self.mass = mass               # kg
        #self.wing_area = wing_area     # m²
        self.chord = chord
        self.wingspan = wingspan       # m
        self.oswald_efficiency = oswald_efficiency
        self.cd0 = cd0                 # drag coefficient at zero lift

    @property
    def wing_area(self):
        return self.chord * (self.wingspan*0.4)
    @property
    def aspect_ratio(self):
        """AR = b² / S"""
        return self.wingspan**2 / self.wing_area

    @property
    def induced_drag_factor(self):
        """K = 1 / (π * AR * e)"""
        return 1 / (math.pi * self.aspect_ratio * self.oswald_efficiency)

    def minimum_sink_velocity(self, rho=1.225):
        """
        V_min_sink ≈ sqrt(2/rho) * (K/CD0)^(1/4) * sqrt(W/S)
        """
        W = self.mass * 9.81  # Newtons
        term_1 = math.sqrt(2 / rho)
        term_2 = (self.induced_drag_factor / self.cd0) ** 0.25
        term_3 = math.sqrt(W / self.wing_area)

        return term_1 * term_2 * term_3
    
glider = Glider(
    mass=0.1,              # 180 g
    chord=0.1,         # m² (ex: 1.2 m x 0.23 m de corda média)
    wingspan=0.8,           # m
    oswald_efficiency=0.75, # típico para asa simples
    cd0=0.028               # arrasto zero-lift de fuselagem+asa básica
)

print("Aspect Ratio:", glider.aspect_ratio)
print("Induced Drag Factor (K):", glider.induced_drag_factor)
print("Velocidade de mínimo sink (m/s):", glider.minimum_sink_velocity())