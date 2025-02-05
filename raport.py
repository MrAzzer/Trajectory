import numpy as np
import matplotlib.pyplot as plt
from numpy._typing import _8Bit

# Stałe
m = 50  # masa w kg
v0 = 100  # prędkość początkowa w m/s
h0 = 100  # wysokość początkowa w m
g = 9.81  # przyspieszenie ziemskie w m/s^2
k = 0.5  # stała oporu powietrza

# Funkcja obliczająca siłę oporu powietrza z uwzględnieniem kierunku prędkości
def f_o(v):
    return k * v * abs(v)

# Metoda Rungego-Kutty
def rungego_kutty(dt=0.0001):
    h, v, t = h0, v0, 0
    h_max = h0  # Inicjalizacja maksymalnej wysokości

    while True:
        k1_v = -g - f_o(v)/m
        k1_h = v
        k2_v = -g - f_o(v + 0.5 * dt * k1_v)/m
        k2_h = v + 0.5 * dt * k1_v
        k3_v = -g - f_o(v + 0.5 * dt * k2_v)/m
        k3_h = v + 0.5 * dt * k2_v
        k4_v = -g - f_o(v + dt * k3_v)/m
        k4_h = v + dt * k3_v

        v_new = v + (dt / 6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        h_new = h + (dt / 6) * (k1_h + 2*k2_h + 2*k3_h + k4_h)

        # Aktualizacja maksymalnej wysokości po obliczeniu nowego h
        h_max = max(h_max, h_new)

        if h_new < 0:
            break

        h, v = h_new, v_new
        t += dt

    return h_max, abs(v_new)

# Metoda Velocity Verlet
def verleta(dt=0.0001):
    h, v, t = h0, v0, 0
    a = -g - f_o(v) / m
    h_max = h0

    while True:
        # Oblicz nowe położenie
        h_new = h + v * dt + 0.5 * a * dt**2

        # Oblicz prędkość w połowie kroku
        v_half = v + 0.5 * a * dt

        # Oblicz nowe przyspieszenie
        a_new = -g - f_o(v_half) / m

        # Oblicz nową prędkość
        v_new = v + 0.5 * (a + a_new) * dt

        # Aktualizacja maksymalnej wysokości
        h_max = max(h_max, h_new)

        if h_new < 0:
            break

        h, v, a = h_new, v_new, a_new
        t += dt

    return h_max, abs(v_new)

# Metoda Heuna
def heuna(dt=0.0001):
    h, v, t = h0, v0, 0
    h_max = h0

    while True:
        a = -g - f_o(v) / m
        v_pred = v + a * dt
        h_pred = h + v * dt
        a_pred = -g - f_o(v_pred) / m

        v_new = v + 0.5 * dt * (a + a_pred)
        h_new = h + 0.5 * dt * (v + v_pred)

        h_max = max(h_max, h_new)

        if h_new < 0:
            break

        h, v = h_new, v_new
        t += dt

    return h_max, abs(v_new)

# Wyniki
max_height_rk, impact_speed_rk = rungego_kutty()
max_height_verlet, impact_speed_verlet = verleta()
max_height_heun, impact_speed_heun = heuna()

print(f"Metoda Rungego-Kutty: Wysokość = {max_height_rk:.8f} m, Prędkość zderzenia = {impact_speed_rk:.8f} m/s")
print(f"Metoda Verleta: Wysokość = {max_height_verlet:.8f} m, Prędkość zderzenia = {impact_speed_verlet:.8f} m/s")
print(f"Metoda Heuna: Wysokość = {max_height_heun:.8f} m, Prędkość zderzenia = {impact_speed_heun:.8f} m/s")
