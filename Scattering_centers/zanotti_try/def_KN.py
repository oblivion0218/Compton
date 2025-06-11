import numpy as np 
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR
import interaction_lib as ilib

file_path = '/mnt/c/Users/User/Desktop/info/Compton/' # Andrea


# DATA
def leggi_dati(file_path):
    # Inizializza liste per ogni colonna
    angle, rate, err_rate = [], [], []
    count, err_count = [], []
    channel, err_channel = [], []
    sigma, err_sigma = [], []

    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Salta l'intestazione
        for line in lines[1:]:
            valori = line.strip().split()
            if len(valori) != 9:
                continue  # Salta righe non valide

            angle.append(float(valori[0]))
            rate.append(float(valori[1]))
            err_rate.append(float(valori[2]))
            count.append(float(valori[3]))
            err_count.append(float(valori[4]))
            channel.append(float(valori[5]))
            err_channel.append(float(valori[6]))
            sigma.append(float(valori[7]))
            err_sigma.append(float(valori[8]))

    return angle, rate, err_rate, count

file_path_reflection = file_path + "Codes/data_analysis/parameters_pol4_riflection.txt"
angle_reflection, rate_reflection, err_rate_reflection, count_reflection = leggi_dati(file_path_reflection)
angle_reflection = np.array(angle_reflection[1:])
rate_reflection = np.array(rate_reflection[1:])
err_rate_reflection = np.array(err_rate_reflection[1:])
count_reflection = np.array(count_reflection[1:])

file_path_trasmission = file_path + "Codes/data_analysis/parameters_pol4_trasmission.txt"
angle_trasmission, rate_trasmission, err_rate_trasmission, count_trasmission = leggi_dati(file_path_trasmission)
angle_trasmission = np.array(angle_trasmission)
rate_trasmission = np.array(rate_trasmission)
err_rate_trasmission = np.array(err_rate_trasmission)
count_trasmission = np.array(count_trasmission)

solid_angles_reflection = np.array([0.02735079041747454, 0.025819652026622272, 0.027070509615309143, 
                                    0.026302777950659273, 0.026575610155933342, 0.02669072072779364,
                                    0.026320830435048827, 0.02677502814759899])
solid_angles_reflection = solid_angles_reflection[1:]
solid_angles_trasmission = np.array([0.026686807818702276, 0.02735079041747454, 
                                     0.025819652026622272, 0.027070509615309143])

# DATA ERRORS
error_reflection = np.load(file_path + 'Measurments/errori_arrays_riflex.npz')
y_err_reflection = np.array(error_reflection['y_err'][1:])
x_err_reflection = np.array(error_reflection['x_err'][1:])

error_trasmission = np.load(file_path + 'Measurments/errori_arrays_trasm.npz')
y_err_trasmission = np.array(error_trasmission['y_err'])
x_err_trasmission = np.array(error_trasmission['x_err'])

statistical_error_reflection = np.array([np.sqrt(err_rate ** 2 + (1/(N))) for err_rate, N in zip(err_rate_reflection, count_reflection)])
statistical_error_trasmission = np.array([np.sqrt(err_rate ** 2 + (1/(N))) for err_rate, N in zip(err_rate_trasmission, count_trasmission)])

y_err_reflection += statistical_error_reflection
y_err_trasmission += statistical_error_trasmission

y_err = np.concatenate((y_err_trasmission, y_err_reflection))
x_err = np.concatenate((x_err_trasmission, x_err_reflection))

# FIT MODELS
def compton_scattering(theta):
    return 511 / (2 - np.cos(np.deg2rad(theta)))  # Compton scattering formula in keV


def efficiency(theta):
    A = 1.5604
    B = -0.0995
    C = 3.53466
    D = 0.10210  
    E = compton_scattering(theta) * 1e-3  # Convert from MeV to keV

    efficiency_value = A * pow(E, -B) * np.exp(-C * E) + D
    return efficiency_value


def constant():
    r_gate = 1.27  #cm
    d_gate_source = 18.54 #cm
    beta_2 = np.arctan(r_gate / d_gate_source)
    Omega =  (1 - np.cos(beta_2)) * 2 * np.pi #sr   

    S = 188900 #Bq
    BR = 0.903
    epsilon_gate = 0.147
    flux = (2 * S * BR) * (Omega / (4 * np.pi)) * epsilon_gate #s^-1

    n_c = ilib.rho * ilib.Z * (ilib.N_a / ilib.MM) 

    return flux * n_c


def KN(parameters, theta):
    A = parameters              # model parameters: amplitude A and initial photon energy E0
    theta = np.deg2rad(theta)        # convert degrees to radians
    E0 = 511.0                       # electron rest energy [keV]
    E_prime = compton_scattering(theta)  # Compton scattering energy in keV
    r = E_prime / E0  # energy ratio

    return A * r**2 * (r + 1/r - np.sin(theta)**2)


def KN_free(parameters, theta):
    A, m_e = parameters              # model parameters: amplitude A and initial photon energy E0
    theta = np.deg2rad(theta)        # convert degrees to radians
    E0 = 511.0                       # electron rest energy [keV]
    # Compute energy ratio E'/E = 1 / [1 + (E0/me)*(1 - cosθ)]
    E_prime = 1.0 / (1 + (E0 / m_e) * (1 - np.cos(theta)))
    r = E_prime / E0  # energy ratio

    return A * (r**2 * (r + 1/r - np.sin(theta)**2))


target_width = 1 # cm
Z = 29 # Atomic number of Copper (Cu)


def prob_trasmission(theta):
    D_prime = target_width
    E = 511
    E_prime = compton_scattering(theta)
    lam = ilib.attenuation_length(E, Z)
    lam_prime = ilib.attenuation_length(E_prime, Z) * np.cos(np.deg2rad(theta))
    lam_dubprime = (lam * lam_prime) / (lam_prime - lam)

    return lam_dubprime * np.exp(- D_prime / lam_prime) * (1 - np.exp(-D_prime / lam_dubprime))


def prob_reflection(theta):
    D_prime = target_width / np.cos((np.pi - np.deg2rad(theta))/2)
    E = 511
    E_prime = compton_scattering(theta)
    lam = ilib.attenuation_length(E, Z)
    lam_prime = ilib.attenuation_length(E_prime, Z)
    lam_dubprime = (lam * lam_prime) / (lam_prime + lam)
    
    return lam_dubprime * (1 - np.exp(-D_prime / lam_dubprime))


# FIT
rate_reflection = [rate_reflection[i]/(prob_reflection(angle) * efficiency(angle)) 
                            for i, angle in enumerate(angle_reflection)]
rate_reflection = rate_reflection / solid_angles_reflection
rate_trasmission = [rate_trasmission[i]/(prob_trasmission(angle) * efficiency(angle))
                             for i, angle in enumerate(angle_trasmission)]
rate_trasmission = rate_trasmission / solid_angles_trasmission
rates = np.concatenate((rate_trasmission, rate_reflection))

angles = np.concatenate((angle_trasmission, angle_reflection))

kn = Model(KN)
data = RealData(angles, rates, sx=x_err, sy=y_err)
odr = ODR(data, kn, beta0=[(constant() * (ilib.r_e ** 2)/2)]) 
out = odr.run()

print("Fitted A =", out.beta[0], "±", out.sd_beta[0])
print("Fitted r_e (cm) =", np.sqrt((2/constant()) * out.beta[0]), "±", np.sqrt((2/constant()) * out.sd_beta[0]))
print("Reduced chi-square =", out.res_var)

kn_free = Model(KN_free)
data_free = RealData(angles, rates, sx=x_err, sy=y_err)
odr_free = ODR(data_free, kn_free, beta0=[(constant() * (ilib.r_e ** 2)/2), ilib.m_e])
out_free = odr_free.run()

print("Fitted A =", out_free.beta[0], "±", out_free.sd_beta[0])
print("Fitted m_e (MeV/c^2) =", out_free.beta[1], "±", out_free.sd_beta[1])
print("Reduced chi-square =", out_free.res_var)

plt.figure(figsize=(10, 6))
plt.errorbar(angles, rates, xerr=x_err, yerr=y_err, fmt='o', label='Data', color='blue')
theta_fit = np.linspace(0, 140, 100)
plt.plot(theta_fit, KN(out.beta, theta_fit), label='Fitted Model', color='red')
plt.plot(theta_fit, KN_free(out_free.beta, theta_fit), label='Fitted Model (Free)', color='green')
plt.xlabel('Angle (degrees)')
plt.ylabel('Rate (normalized)')
plt.title('Compton Scattering Fit')
plt.legend()
plt.grid()
plt.savefig('Compton_scattering_fit.png')

