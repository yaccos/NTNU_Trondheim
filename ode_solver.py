import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(y, t, k):
    dydt = - y * k
    return dydt

def test():
    y0 = 5
    t = np.linspace(0,20)
    k = 0.1
    y = odeint(model, y0, t, args=(k,))
    
    plot_tool.plot_ode(t,y)

def test_HIV():
    # https://apmonitor.com/pdc/index.php/Main/SimulateHIV
    t = np.linspace(0,15,5*365)
    healthy = H0 = 1e6
    infected = I0 = 0
    virus = V0 = 100
    x0 = [H0, I0, V0]

    x = odeint(model_HIV, x0, t)
    
    H = x[:,0]; I = x[:,1]; V = x[:,2]
    plt.semilogy(t,H,label="Healthy")
    plt.semilogy(t,I,label="Infected")
    plt.semilogy(t,V,label="Virus")
    plt.legend()
    plt.show()


def model_HIV(x,t):
    # https://apmonitor.com/pdc/index.php/Main/SimulateHIV
    H, I, V = x[0], x[1], x[2]
    k1, k2, k3, k4, k5, k6 = 1e5, 0.1, 2e-7, 0.5, 5, 100
    dHdt = k1 - k2 * H - k3*H*V
    dIdt = k3 * H * V - k4 * I
    dVdt = -k3 * H * V - k5 * V + k6 * I
    return [dHdt, dIdt, dVdt]

test_HIV()
