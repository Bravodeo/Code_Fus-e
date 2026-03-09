import math
import matplotlib.pyplot as plt

def simulate_water_rocket(

    volume_bottle=1.5e-3,
    mass_dry = 0.256, # Masse fusée
    mass_water_init = 0.5, # Masse d'eau initial
    p_air_init = 5.0 * 101325.0, # Pression initial
    p_atm = 101325.0, # Pression athmosphèrique
    gamma = 1.4, #Constante adiabatique
    Ae = 1e-4, #Air tuyère
    rho_water = 1000.0, #Masse volumique de l'eau
    Cd = 0.5, #Coefficient de traînée
    A_cross = 0.01, # Air en contact avec l'air
    rho_air = 1.225, #Masse volumique de l'air
    v0 = 0.0, #vitesse intial
    h0 = 0.0, #Hauteur initial
    dt = 1e-4, #Pas de temps
    max_time = 60.0 #Temps maximal étudié


):
    g = 9.80665 # Constante gravitationel
    t = 0.0 #Temp initial
    v = float(v0)
    h = float(h0)
    m_w = float(mass_water_init)
    m = mass_dry + m_w # Masse total

    V_air_init = volume_bottle - (mass_water_init / rho_water)
    if V_air_init <= 0:
        raise ValueError("Volume d'air initial négatif ou nul.")

    const_pV = p_air_init * (V_air_init ** gamma)

    times = []
    heights = []
    speeds = []
    thrusts = []
    pressures = []
    masses = []

    thrust_phase = True

    while t < max_time:
        V_air = volume_bottle - (m_w / rho_water)
        p_air = const_pV / (V_air ** gamma)

        if thrust_phase and (m_w > 1e-8) and (p_air > p_atm):
            delta_p = p_air - p_atm
            if delta_p <= 0:
                v_e = 0.0
            else:
                v_e = math.sqrt(2.0 * delta_p / rho_water)

            dm_w_dt = - rho_water * Ae * v_e
            dm_w = dm_w_dt * dt

            if (m_w + dm_w) < 0:
                dm_w = -m_w
                dm_w_dt = dm_w / dt

            thrust = rho_water * Ae * (v_e ** 2)
            m_w = m_w + dm_w
            m = mass_dry + m_w
        else:
            thrust = 0.0
            v_e = 0.0
            dm_w_dt = 0.0
            thrust_phase = False

        drag = 0.5 * rho_air * Cd * A_cross * v * abs(v)
        weight = m * g
        a = (thrust - drag - weight) / m

        v = v + a * dt
        h = h + v * dt
        t = t + dt

        times.append(t)
        heights.append(h)
        speeds.append(v)
        thrusts.append(thrust)
        pressures.append(p_air)
        masses.append(m)

        if (not thrust_phase) and (v <= 0):
            break

    # Détection fin de poussée
    t_burn = None
    h_burn = None
    v_burn = None

    for i in range(len(times)):
        if thrusts[i] == 0.0:
            idx = max(0, i - 1)
            t_burn = times[idx]
            h_burn = heights[idx]
            v_burn = speeds[idx]
            break

    if t_burn is None:
        t_burn = 0.0
        h_burn = heights[0] if len(heights) > 0 else h0
        v_burn = speeds[0] if len(speeds) > 0 else v0

    return {
        "time": times,
        "height": heights,
        "speed": speeds,
        "thrust": thrusts,
        "pressure": pressures,
        "mass": masses,
        "t_burn_end": t_burn,
        "h_burn_end": h_burn,
        "v_burn_end": v_burn,
        "t_apogee": t,
        "h_apogee": h,
        "v_apogee": v
    }




res = simulate_water_rocket()
max_thrust = max(res["thrust"])
print("Poussée maximale :", max_thrust, "N")
idx = res["thrust"].index(max_thrust)
t_max_thrust = res["time"][idx]


print("Temps fin poussée :", res["t_burn_end"])
print("Hauteur fin poussée :", res["h_burn_end"])
print("Vitesse fin poussée :", res["v_burn_end"])
print("Hauteur max :", res["h_apogee"])
print("Temps apogée :", res["t_apogee"])


plt.figure(figsize=(8,5))
plt.plot(res["time"], res["height"])
plt.xlabel("Temps (s)")
plt.ylabel("Hauteur (m)")
plt.title("Hauteur en fonction du temps - Fusée à eau")
plt.grid(True)
plt.show()

