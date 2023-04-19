import numpy as np

def model_1(X, p, b):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / ((2 * xi * (1 + gamma) ** 2) * \
        (ie - ea)) - 2494.35) / (v_m)) + b
    
def model_2(X, p, b):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((ie + ea) ** 2 / (8 * xi * (ie - ea)) - 2494.35) / \
        (v_m)) + b
    
def model_3(X, p, b):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / ((2 * (1 + gamma) ** 2) * \
        (ie - ea)) - 2494.35) / (v_m)) + b

def model_4(X, p, b, h):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((ie + ea) ** 2 / (8 * (ie - ea)) * (1 + ((ie + ea) / \
        (6 * (ie - ea) ** 3) * h)) - 2494.35) / (v_m)) + b

def model_5(X, p, b, a):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((ie + ea) ** 2 / (8 * (ie - ea)) * (1 + ((ie + ea) /\
        (6 * (ie - ea) ** 2) * a)) - 2494.35) / (v_m)) + b

def model_6(X, p, b, h):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / (2 * xi * (1 + gamma) ** 2 * \
        (ie - ea)) * (1 + ((gamma * ie + ea) / (3 * (1 + gamma) * xi ** 3 * \
        (ie - ea) ** 3) * h)) - 2494.35) / (v_m)) + b

def model_7(X, p, b, a):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / (2 * xi * (1 + gamma) ** 2 * \
        (ie - ea)) * (1 + ((gamma * ie + ea) / (3 * (1 + gamma) * xi ** 2 * \
        (ie - ea) ** 2) * a)) - 2494.35) / (v_m)) + b

def model_8(X, p, b, h):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((ie + ea) ** 2 / (8 * xi * (ie - ea)) * (1 + ((ie + ea) / \
        (6 * xi ** 3 * (ie - ea) ** 3) * h)) - 2494.35) / (v_m)) + b

def model_9(X, p, b, a):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / (8 * xi * (ie - ea)) * (1 + ((ie + ea) / \
        (6 * xi ** 2 * (ie - ea) ** 2) * a)) - 2494.35) / (v_m)) + b
    
def model_10(X, p, b, h):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / (2 * (1 + gamma) ** 2 * \
        (ie - ea)) * (1 + ((gamma * ie + ea) / (3 * (1 + gamma) * \
        (ie - ea) ** 3) * h)) - 2494.35) / (v_m)) + b

def model_11(X, p, b, a):
    homo, lumo, v_m, gamma, xi = X
    ie = -homo
    ea = -lumo
    return p * np.sqrt(((gamma * ie + ea) ** 2 / (2 * (1 + gamma) ** 2 * \
        (ie - ea)) * (1 + ((gamma * ie + ea) / (3 * (1 + gamma) * \
        (ie - ea) ** 2) * a)) - 2494.35) / (v_m)) + b
