import math

class Glider:
    def __init__(self, mass, chord, wingspan, oswald_efficiency, cd0):
        self.mass = mass               # kg
        self.chord = chord             # m (corda média)
        self.wingspan = wingspan       # m (envergadura ponta-a-ponta)
        self.oswald_efficiency = oswald_efficiency
        self.cd0 = cd0                 # drag coefficient at zero lift

    @property
    def wing_area(self):
        # área da asa (aproximação rectangular ou média de corda)
        return self.chord * self.wingspan

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
        retorna velocidade de avanço (m/s) que minimiza o sink
        """
        W = self.mass * 9.81  # Newtons
        term_1 = math.sqrt(2 / rho)
        term_2 = (self.induced_drag_factor / self.cd0) ** 0.25
        term_3 = math.sqrt(W / self.wing_area)
        return term_1 * term_2 * term_3

    def vertical_sink_rate(self, V=None, rho=1.225):
        """
        Calcula a taxa de sink vertical (m/s).
        Se V for None usa a velocidade de mínimo sink.
        Também retorna um dicionário com CL, CD, D, glide_ratio para debug.
        """
        W = self.mass * 9.81
        S = self.wing_area

        if V is None:
            V = self.minimum_sink_velocity(rho=rho)

        # coeficiente de sustentação requerido
        CL = W / (0.5 * rho * V**2 * S)

        # coeficiente de arrasto (parcial: CD0 + K*CL^2)
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL**2

        # arrasto (N)
        D = 0.5 * rho * V**2 * S * CD

        # sink vertical (m/s) ~ V * (D / W)
        sink_rate = V * (D / W)

        # glide ratio L/D ≈ W/D (assumindo L~=W em descida estabilizada)
        glide_ratio = W / D if D != 0 else float('inf')

        return {
            "V_used": V,
            "CL": CL,
            "CD": CD,
            "D_N": D,
            "sink_rate_m_s": sink_rate,
            "glide_ratio": glide_ratio
        }

    def stall_velocity(self, CL_max, rho=1.225):
        W = self.mass * 9.81
        
        S = self.wing_area

        return ((2*W)/(rho*S*CL_max))
    
# --- exemplo com os teus valores ---
glider = Glider(
    mass=0.3,              # kg
    chord=0.1,             # m
    wingspan=1,           # m
    oswald_efficiency=0.75, # típico para asa simples
    cd0=0.04                # arrasto zero-lift de fuselagem+asa básica
)

print("Aspect Ratio:", glider.aspect_ratio)
print("Induced Drag Factor (K):", glider.induced_drag_factor)
Vmin = glider.minimum_sink_velocity()
print("Velocidade de mínimo sink (m/s):", Vmin)

sink_info = glider.vertical_sink_rate()  # usa V_min por omissão
print("Usado V (m/s):", sink_info["V_used"])
print("CL:", sink_info["CL"])
print("CD:", sink_info["CD"])
print("Arrasto D (N):", sink_info["D_N"])
print("Taxa de sink vertical (m/s):", sink_info["sink_rate_m_s"])
print("Glide ratio (L/D):", sink_info["glide_ratio"])

print(f'Stall Velocity: {glider.stall_velocity(1.4)} m/s')
