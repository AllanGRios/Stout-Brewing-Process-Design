import matplotlib.pyplot as plt

# Constants
thermal_conductivity = 16.2  # W/(m K)
Cp_ss = 500  # J/(kg K)
Cp_water = 4184  # J/(kg K)
density_ss = 7930  # kg/m3
density_water = 1000  # kg/m3
emissivity_ss304 = 0.8
emissivity_water = 0.963
heat_transfer_coef = 25  # W/(m2 k)
T_env_C = 23  # Celcius
T_env_K = T_env_C + 273.15  # Kelvin
stephan_boltzman_coeff = 5.67e-8  # W/(m2 K4)
A_side = 0.0806  # m2
A_top = 0.046  # m2

T_start_C = 72.0
simulation_minutes = 90
dt = 1
total_steps = simulation_minutes * 60

# total heat capacity
mass_ss = (21.324 / 1e6) * density_ss  # kg (pot mass)
mass_water = 0.0025 * density_water  # kg (2.5l water)
total_heat_capacity = (mass_ss * Cp_ss) + (mass_water * Cp_water)

#plotting lists
time_axis = []
temperature_history = []
fc_loss_history = []
rad_loss_history = []

current_T_C = T_start_C

for step in range(total_steps):
    current_time_seconds = step * dt
    current_T_K = current_T_C + 273.15
    # free convection
    dT = current_T_C - T_env_C
    q_fc = heat_transfer_coef * (dT * A_side + dT * A_top)  # j/s

    # radiation
    q_rad = stephan_boltzman_coeff * ((current_T_K**4 - T_env_K**4) * emissivity_ss304 * A_side + (current_T_K**4 - T_env_K**4) * emissivity_water * A_top)  # j/s

    # Total loss
    total_power_lost = q_fc + q_rad

    # Instantaneous temperature drop
    # Q = P * t   dT = Q / Cp_total
    energy_lost = total_power_lost * dt
    temperature_drop = energy_lost / total_heat_capacity

    current_T_C -= temperature_drop

    time_axis.append(current_time_seconds / 60.0)
    temperature_history.append(current_T_C)
    fc_loss_history.append(q_fc)
    rad_loss_history.append(q_rad)
    if 61.99 < temperature_history[step] < 62:
        SEC = (step/60 - step//60)*60
        MIN = step/60 - (step/60 - step//60)
    elif 61.98 < temperature_history[step] < 62:
        SEC = (step / 60 - step // 60) * 60
        MIN = step / 60 - (step / 60 - step // 60)
plt.figure(figsize=(12, 5))

# Temperature Decay Curve
plt.subplot(1, 2, 1)
plt.plot(time_axis, temperature_history, color="red", label="Uninsulated Pot")
plt.axhline(62, color="blue", linestyle="--", label=f"Minimum Wort Temp (62°C)\nTime to minimum {int(MIN)} min {int(SEC)} s")
plt.title("Wort Temperature Profile over 90 mins")
plt.xlabel("Time (Minutes)")
plt.ylabel("Temperature (°C)")
plt.grid(True)
plt.legend()

# Convection vs Radiation graph
plt.subplot(1, 2, 2)
plt.plot(time_axis, fc_loss_history, label="Convection Loss ( J/s )", color="orange")
plt.plot(time_axis, rad_loss_history, label="Radiation Loss ( J/s )", color="purple")
plt.title("Heat Transfer Mechanism Comparison")
plt.xlabel("Time (Minutes)")
plt.ylabel("Energy Loss Rate ( J/s )")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()