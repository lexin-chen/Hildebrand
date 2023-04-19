import numpy as np

class Models:
    def __init__(self, homo, lumo, v_m, **kwargs):
        self.ie = -homo
        self.ea = -lumo
        self.v_m = v_m
        self.gamma = kwargs.get("gamma", None)
        self.xi = kwargs.get("xi", None)

    def model_1(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / ((2 * self.xi * (1 + self.gamma) ** 2) * \
            (self.ie - ea)) - 2494.35) / (self.v_m)) + b
        
    def model_2(self):
        return p * np.sqrt(((self.ie + ea) ** 2 / (8 * self.xi * (self.ie - ea)) - 2494.35) / \
            (self.v_m)) + b
        
    def model_3(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / ((2 * (1 + self.gamma) ** 2) * \
            (self.ie - ea)) - 2494.35) / (self.v_m)) + b

    def model_4(self):
        return p * np.sqrt(((self.ie + ea) ** 2 / (8 * (self.ie - ea)) * (1 + ((self.ie + ea) / \
            (6 * (self.ie - ea) ** 3) * h)) - 2494.35) / (self.v_m)) + b

    def model_5(self):
        return p * np.sqrt(((self.ie + ea) ** 2 / (8 * (self.ie - ea)) * (1 + ((self.ie + ea) /\
            (6 * (self.ie - ea) ** 2) * a)) - 2494.35) / (self.v_m)) + b

    def model_6(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / (2 * self.xi * (1 + self.gamma) ** 2 * \
            (self.ie - ea)) * (1 + ((self.gamma * self.ie + ea) / (3 * (1 + self.gamma) * self.xi ** 3 * \
            (self.ie - ea) ** 3) * h)) - 2494.35) / (self.v_m)) + b

    def model_7(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / (2 * self.xi * (1 + self.gamma) ** 2 * \
            (self.ie - ea)) * (1 + ((self.gamma * self.ie + ea) / (3 * (1 + self.gamma) * self.xi ** 2 * \
            (self.ie - ea) ** 2) * a)) - 2494.35) / (self.v_m)) + b

    def model_8(self):
        return p * np.sqrt(((self.ie + ea) ** 2 / (8 * self.xi * (self.ie - ea)) * (1 + ((self.ie + ea) / \
            (6 * self.xi ** 3 * (self.ie - ea) ** 3) * h)) - 2494.35) / (self.v_m)) + b

    def model_9(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / (8 * self.xi * (self.ie - ea)) * (1 + ((self.ie + ea) / \
            (6 * self.xi ** 2 * (self.ie - ea) ** 2) * a)) - 2494.35) / (self.v_m)) + b
        
    def model_10(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / (2 * (1 + self.gamma) ** 2 * \
            (self.ie - ea)) * (1 + ((self.gamma * self.ie + ea) / (3 * (1 + self.gamma) * \
            (self.ie - ea) ** 3) * h)) - 2494.35) / (self.v_m)) + b

    def model_11(self):
        return p * np.sqrt(((self.gamma * self.ie + ea) ** 2 / (2 * (1 + self.gamma) ** 2 * \
            (self.ie - ea)) * (1 + ((self.gamma * self.ie + ea) / (3 * (1 + self.gamma) * \
            (self.ie - ea) ** 2) * a)) - 2494.35) / (self.v_m)) + b
       