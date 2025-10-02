# Criar um planador
glider = Glider(
    mass=400,             # kg
    wing_area=16,         # m²
    wingspan=20,          # m
    oswald_efficiency=0.85,
    cd0=0.02
)

print("Aspect Ratio:", glider.aspect_ratio)
print("Induced Drag Factor (K):", glider.induced_drag_factor)
print("Velocidade de mínimo sink (m/s):", glider.minimum_sink_velocity())
