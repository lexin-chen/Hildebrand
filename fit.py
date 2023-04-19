import pandas as pd
import numpy as np
from enum import Enum
from scipy.optimize import curve_fit
import modules.models as models
from matplotlib import pyplot as plt

class prop_col(Enum):
    """
    Enum class for the column number of the properties
    """
    HSP_exp = 1
    HOMO = 2
    LUMO = 3
    V_m = 4
    Dipole = 5
    Hydrogen = 6
    chi = 7
    eta = 8

def get_property_col(property):
    """ Get the column number of the property from the Enum class prop_col
    Args:
        property (str): List of all possible properties.
    Returns:
        int: Column number of the property
    """
    return prop_col[property].value

def read_data(file_name):
    # property = get_property_col(property)
    df = pd.read_csv(file_name, delimiter=",", usecols=[0, 1, 2, 3, 4], 
                            index_col=0, encoding='utf-8')

if __name__ == "__main__":
    file_name = "all_molecules.csv"
    df = pd.read_csv(file_name, delimiter=",", encoding='utf-8')
    func = getattr(models, "model_7")
    parameters, covariance = curve_fit(func, (df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), df['HSP_exp'], maxfev = 500000)
    fit_data = func((df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), *parameters)
    print(fit_data)
    SE = np.sqrt(np.diag(covariance))
    fig, ax = plt.subplots()
    x = np.linspace(min(df['HSP_exp']), max(df['HSP_exp']))
    ax.scatter(df['HSP_exp'], fit_data, marker='o')
    ax.plot(x, x, dashes=[4, 4])
    plt.show()