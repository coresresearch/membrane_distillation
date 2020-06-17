"""
Defines functions for running membrane distillation simulation.

- Residual_1d: Calcualtes time-based derivative of solution vector, dSV/dt, for 
    a 1-D simulation.
"""
import numpy as np

def residual_1d(t, SV, obj, params):
    residual = np.zeros_like(SV)

    return residual