import math  # Bibliothèque mathématique (racine carrée, etc.)
import matplotlib.pyplot as plt  # Bibliothèque pour tracer des graphiques

# Fonction principale qui simule le vol d'une fusée à eau
def simulate_water_rocket(

    volume_bottle=1.5e-3,  # Volume total de la bouteille (m³)
    mass_dry = 0.256, # Masse de la fusée vide (kg)
    mass_water_init = 0.6, # Masse initiale d'eau dans la fusée (kg)
    p_air_init = 5.0 * 100000.0, # Pression initiale de l'air dans la fusée (Pa)
    p_atm = 101325.0, # Pression atmosphérique extérieure (Pa)
    gamma = 1.4, # Constante adiabatique de l'air (rapport Cp/Cv)
    Ae = 9e-5, # Surface de la tuyère (m²)
    rho_water = 1000.0, # Masse volumique de l'eau (kg/m³)
    Cd = 0.5, # Coefficient de traînée aérodynamique
    A_cross = 0.01, # Surface frontale de la fusée en contact avec l'air (m²)
    rho_air = 1.225, # Masse volumique de l'air (kg/m³)
    C_discharge = 0.8, # Coefficient de décharge de la tuyère (efficacité de l'écoulement)
    v0 = 0.0, # Vitesse initiale de la fusée (m/s)
    h0 = 0.0, # Hauteur initiale de la fusée (m)
    dt = 1e-4, # Pas de temps de la simulation (s)
    max_time = 60.0 # Temps maximum simulé (s)

):
    g = 9.80665 # Accélération gravitationnelle terrestre (m/s²)
    t = 0.0 # Temps initial
    v = float(v0) # Vitesse initiale convertie en float
    h = float(h0) # Hauteur initiale convertie en float
    m_w = float(mass_water_init) # Masse actuelle d'eau
    m = mass_dry + m_w # Masse totale de la fusée (structure + eau)

    # Calcul du volume d'air initial dans la bouteille
    V_air_init = volume_bottle - (mass_water_init / rho_water)

    # Vérification : il faut qu'il reste de l'air dans la bouteille
    if V_air_init <= 0:
        raise ValueError("Volume d'air initial négatif ou nul.")

    # Constante de la loi adiabatique pV^gamma = constante
    const_pV = p_air_init * (V_air_init ** gamma)

    # Listes servant à stocker les résultats au cours du temps
    times = []
    heights = []
    speeds = []
    thrusts = []
    pressures = []
    masses = []

    thrust_phase = True # Indique si la phase de poussée est encore active

    # Boucle principale de simulation temporelle
    while t < max_time:

        # Calcul du volume d'air actuel dans la bouteille
        V_air = volume_bottle - (m_w / rho_water)

        # Calcul de la pression de l'air via la loi adiabatique
        p_air = const_pV / (V_air ** gamma)

        # Vérifie si la fusée produit encore de la poussée
        if thrust_phase and (m_w > 1e-8) and (p_air > p_atm):

            # Différence de pression entre l'intérieur et l'extérieur
            delta_p = p_air - p_atm

            # Si la pression est insuffisante, pas d'éjection
            if delta_p <= 0:
                v_e = 0.0
            else:
                # Vitesse d'éjection de l'eau (équation de Bernoulli)
                v_e = math.sqrt(2.0 * delta_p / rho_water)

            # Débit massique d'eau sortant par la tuyère
            dm_w_dt = - rho_water * (C_discharge * Ae) * v_e

            # Masse d'eau éjectée pendant le pas de temps
            dm_w = dm_w_dt * dt

            # Empêche que la masse d'eau devienne négative
            if (m_w + dm_w) < 0:
                dm_w = -m_w
                dm_w_dt = dm_w / dt

            # Calcul de la poussée générée par l'éjection de l'eau
            thrust = rho_water * (C_discharge * Ae) * (v_e ** 2)

            # Mise à jour de la masse d'eau restante
            m_w = m_w + dm_w

            # Mise à jour de la masse totale
            m = mass_dry + m_w

        else:
            # Si plus de poussée
            thrust = 0.0
            v_e = 0.0
            dm_w_dt = 0.0
            thrust_phase = False

        # Calcul de la force de traînée aérodynamique
        drag = 0.5 * rho_air * Cd * A_cross * v * abs(v)

        # Force du poids
        weight = m * g

        # Accélération résultante (Seconde loi de Newton)
        a = (thrust - drag - weight) / m

        # Mise à jour de la vitesse
        v = v + a * dt

        # Mise à jour de la position (hauteur)
        h = h + v * dt

        # Mise à jour du temps
        t = t + dt

        # Sauvegarde des données pour analyse
        times.append(t)
        heights.append(h)
        speeds.append(v)
        thrusts.append(thrust)
        pressures.append(p_air)
        masses.append(m)

        # Arrêt de la simulation lorsque la fusée commence à redescendre
        if (not thrust_phase) and (v <= 0):
            break

    # Détection de la fin de la poussée
    t_burn = None
    h_burn = None
    v_burn = None

    # Recherche du moment où la poussée devient nulle
    for i in range(len(times)):
        if thrusts[i] == 0.0:
            idx = max(0, i - 1)
            t_burn = times[idx]
            h_burn = heights[idx]
            v_burn = speeds[idx]
            break

    # Cas où la poussée n'est pas détectée
    if t_burn is None:
        t_burn = 0.0
        h_burn = heights[0] if len(heights) > 0 else h0
        v_burn = speeds[0] if len(speeds) > 0 else v0

    # Retourne tous les résultats de la simulation dans un dictionnaire
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


# Exécution de la simulation avec les paramètres par défaut
res = simulate_water_rocket()

# Recherche de la poussée maximale
max_thrust = max(res["thrust"])
print("Poussée maximale :", max_thrust, "N")

# Trouve l'instant correspondant à cette poussée maximale
idx = res["thrust"].index(max_thrust)
t_max_thrust = res["time"][idx]

# Affiche différentes informations sur le vol
print("Temps fin poussée :", res["t_burn_end"])
print("Hauteur fin poussée :", res["h_burn_end"])
print("Vitesse fin poussée :", res["v_burn_end"])
print("Hauteur max :", res["h_apogee"])
print("Temps apogée :", res["t_apogee"])


# Création d'une figure pour afficher le graphique
plt.figure(figsize=(8,5))

# Trace la hauteur en fonction du temps
plt.plot(res["time"], res["height"])

# Label de l'axe des x
plt.xlabel("Temps (s)")

# Label de l'axe des y
plt.ylabel("Hauteur (m)")

# Titre du graphique
plt.title("Hauteur en fonction du temps - Fusée à eau")

# Affiche une grille
plt.grid(True)

# Affiche le graphique
plt.show()
