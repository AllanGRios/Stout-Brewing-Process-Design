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

batch1 = Batch("Batch 1", 1.044, 1.010, 27, 1.7185, 0.45, 0.15, 0.085, 0.11, 1.164, 1.5, 3)
batch2 = Batch("Batch 2", 1.042, 1.008, 27, 2.729, 0.45, 0.15, 0.085, 0.11, 0.73, 2.2, 3)

def data_processing(batch):
    OG = batch.OG
    FG = batch.FG
    C0 = 2714 * OG - 2713
    C1 = 2714 * FG - 2713
    P = C0 - C1
    return C0, P
def Gompertz_model(time, batch):
    P_t = P * np.exp( -1 * np.exp(((batch.Rmax * np.e)/(P)) * (batch.lag - time) + 1))
    return P_t

def Monod_Growth_model(time, batch):
    C = C0
    X = batch.pitching_rate

    C_history = []
    X_histroy = []
    mu_history = []

    dt = 1

    for t in time:
        mu = batch.mu_max * ((C)/(batch.Ks + C)) * np.exp(-1 * batch.Ki * 0.51143 * (C0 - C))
        dX = mu * X * dt
        dC = (- 1/batch.Yxs) * dX

        C += dC
        X += dX

        C_history.append(C0 - C)
        X_histroy.append(X)
        mu_history.append(mu)

    return C_history, X_histroy, mu_history

def main(batch):
    timescale = 72
    time = np.linspace(0,timescale, timescale + 1)
    P_t = []
    for i in time:
        y1 = Gompertz_model(i, batch)
        P_t.append(y1)

    sugar, biomass, mu = Monod_Growth_model(time, batch)
    
    plt.figure(figsize=(10,7))

    plt.subplot(2, 1, 1)
    plt.plot(time,P_t, label=f"{batch.name} \n Sugar Consumed {round(((max(P_t) - min(P_t)) * batch.V ), 1)} g")
    plt.xlim(min(time))
    plt.ylim(min(P_t))
    plt.title("Modified Gompertz For Sugar Consumption")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Sugar Consumed (g/L)")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(time, sugar, label=f"{batch.name} \n Sugar Consumed {round(((max(sugar) - min(sugar)) * batch.V ), 1)} g")
    plt.xlim(min(time))
    plt.ylim(min(sugar))
    plt.title("Monod Kinetic Model")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Sugar Consumed (g/L)")
    plt.grid(True)
    plt.legend()

    plt.subplots_adjust(wspace=0.1, hspace=0.5)
    plt.show()

command = int(input(f"Which batch are you looking for? \nType in the corresponding batch number - "))
corresponding = {1: batch1, 2: batch2}

C0, P = data_processing(corresponding[command])
main(corresponding[command])