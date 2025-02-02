import numpy as np
import matplotlib.pyplot as plt

# Stałe
m = 50  # masa w kg
v0 = 100  # prędkość początkowa w m/s
h0 = 100  # wysokość początkowa w m
g = 9.81  # przyspieszenie ziemskie w m/s^2
k = 0.5  # stała oporu powietrza

# Funkcja obliczająca siłę oporu powietrza
def f_o(v):
    return k * v*v

# Metoda Rungego-Kutty
def rungego_kutty(dt=0.0001):
    h, v, t = h0, v0, 0 #definicja wysokosci, predkosci i czasu
    #Dopóki obiekt nie dotknie ziemi
    while h >= 0:
        #Obliczanie wartości pochodnej funkcji v(t) oraz h(t) dla kroków k1, k2, k3, k4
        k1_v = -g - f_o(v)/m
        k1_h = v

        k2_v = -g - f_o(v + 0.5 * dt * k1_v)/m
        k2_h = v + 0.5 * dt * k1_v

        k3_v = -g - f_o(v + 0.5 * dt * k2_v)/m
        k3_h = v + 0.5 * dt * k2_v

        k4_v = -g - f_o(v + dt * k3_v)/m
        k4_h = v + dt * k3_v

        #Obliczanie nowych wartości funkcji v(t) oraz h(t)
        v += (dt / 6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        h += (dt / 6) * (k1_h + 2*k2_h + 2*k3_h + k4_h)
        t += dt #aktualizacja wartości kroku czasu

    #Zwracanie maksymalnej wysokości i prędkości przy dotknięciu ziemi
    return h0 + (v0**2) / (2 * g), abs(v)


# Metoda Verleta
def verleta(dt=0.0001):
    h, v, t = h0, v0, 0 #definicja wysokosci, predkosci i czasu
    a = -g - f_o(v) / m #obliczenie przyspieszenia
    h_prev = h - v * dt + 0.5 * a * dt**2 #obliczenie poprzedniej wysokosci
    while h >= 0:
        #Obliczanie nowej pozycji
        h_new = 2 * h - h_prev + a * dt**2
        #Obliczanie nowej prędkości
        v = (h_new - h_prev) / (2 * dt)
        #obliczanie nowego przyspieszenia
        a = -g - f_o(v) / m
        #aktualizacja kroku i powrócenie do początku pętli
        h_prev, h = h, h_new
        t += dt

    return h0 + (v0**2) / (2 * g), abs(v)

# Metoda Heuna
def heuna(dt=0.0001):
    # Inicjalizacja zmiennych wysokosci, prędkości i czasu
    h, v, t = h0, v0, 0

    while h >= 0:
        # Obliczanie przyspieszenia
        a = -g - f_o(v) / m
        # Obliczanie przewidywanej pozycji
        v_pred = v + a * dt
        # Obliczanie przewidywanej prędkości
        h_pred = h + v * dt
        # Obliczanie przewidywanego przyspieszenia
        a_pred = -g - f_o(v_pred) / m
        # iteracyjne obliczanie nowych wartości funkcji v(t) oraz h(t)
        v += 0.5 * dt * (a + a_pred)
        h += 0.5 * dt * (v + v_pred)
        t += dt

    return h0 + (v0**2) / (2 * g), abs(v)

# Wyniki
max_height_rk, impact_speed_rk = rungego_kutty()
max_height_verlet, impact_speed_verlet = verleta()
max_height_heun, impact_speed_heun = heuna()

print(f"Metoda Rungego-Kutty: Wysokość = {max_height_rk:.8f} m, Prędość zderzenia = {impact_speed_rk:.8f} m/s")
print(f"Metoda Verleta: Wysokość = {max_height_verlet:.8f} m, Prędość zderzenia = {impact_speed_verlet:.8f} m/s")
print(f"Metoda Heuna: Wysokość = {max_height_heun:.8f} m, Prędość zderzenia = {impact_speed_heun:.8f} m/s")

