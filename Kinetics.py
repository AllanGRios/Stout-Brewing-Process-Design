import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass

@dataclass
class Batch:
    name: str
    OG: float
    FG: float
    temperature: float # C Max Fermentation
    V: float # L Batch Volume
    mu_max: float # 1/h Max growth rate Voss Kveik
    Ks: float # g/L Half sat const
    Ki: float # L/g Ethanol Inhibition Const
    Yxs: float # g/g Biomass Yield Coeff
    pitching_rate: float # g/L
    lag: float # hrs apparent lag time
    Rmax: float # g/(L h) Max rate of sugar consumption

batch1 = Batch("Batch 1", 1.044, 1.010, 27, 1.7185, 0.45, 0.025, 0.015, 0.1, 1.164, 1.5, 2.8494)
batch2 = Batch("Batch 2", 1.042, 1.008, 27, 2.729, 0.45, 0.025, 0.015, 0.1, 0.73, 2.2, 1.1417)

def data_processing(batch):
    OG = batch.OG
    FG = batch.FG
    C0 = 2714 * OG - 2713
    C1 = 2714 * FG - 2713
    P = C0 - C1
    ABV = ( OG - FG ) * 131.25
    E1 = 7.89 * ABV
    Yes = E1 / P
    return C0, P, E1, Yes
def Gompertz_model(time, batch):
    E_max = E1
    S_t = S_max * np.exp(-1 * np.exp(((batch.Rmax * np.e) / (S_max)) * (batch.lag - time) + 1))
    E_t = E_max * np.exp(-1 * np.exp(((batch.Rmax * np.e) / (E_max)) * (batch.lag - time) + 1))
    return S_t, E_t
def Monod_Growth_model(time, batch):
    C = C0
    X = batch.pitching_rate
    E = 0

    C_history = []
    X_history = []
    E_history = []

    dt = time[1] - time[0]

    for t in time:
        I = (1 - E / 94.71) ** 2
        mu = batch.mu_max * (C/(batch.Ks + C)) * I # Yeast growth rate
        dX = mu * X * (1 - X/8) * dt
        R = batch.Rmax * C / (batch.Ks + C) * I
        dC = - R * X * dt
        dE = - Yes * dC

        C += dC
        E += dE
        X += dX
        C = max(C, 0)
        C_history.append(C0-C)
        E_history.append(E)
        X_history.append(X)
    return C_history, X_history, E_history
def main(batch, timescale=24):
    time = np.linspace(0,timescale, timescale + 1)
    S_t = []
    E_t = []
    for i in time:
        y1, y2 = Gompertz_model(i, batch)
        S_t.append(y1)
        E_t.append(y2)

    sugar, biomass, E = Monod_Growth_model(time, batch)
    
    plt.figure(figsize=(10,7))

    plt.subplot(2, 1, 1)
    plt.plot(time,S_t, label=f"{batch.name} \n Sugar Consumed {round(((max(S_t) - min(S_t)) * batch.V ), 1)} g")
    plt.plot(time, E_t, label=f"{batch.name} Ethanol Produced")
    plt.xlim(min(time))
    plt.ylim(min(S_t))
    plt.title("Modified Gompertz Model")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Concentration (g/L)")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(time, sugar, label=f"{batch.name} \n Sugar Consumed {round(((max(sugar) - min(sugar)) * batch.V ), 1)} g")
    plt.plot(time, E, label=f"{batch.name} Ethanol Produced")
    plt.xlim(min(time))
    plt.ylim(min(E))
    plt.title("Modified Monod-Logistic Model")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Concentration (g/L)")
    plt.grid(True)
    plt.legend()

    plt.subplots_adjust(wspace=0.1, hspace=0.5)
    plt.show()

command = int(input(f"Which batch are you looking for? \nType in the corresponding batch number - "))
corresponding = {1: batch1, 2: batch2}
requested_timescale = int(input(f"What modelling timescale are you looking for? \n Type a number in hours - "))

C0, S_max, E1, Yes = data_processing(corresponding[command])
main(corresponding[command], requested_timescale)