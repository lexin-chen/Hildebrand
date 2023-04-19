import pandas as pd
import numpy as np
from enum import Enum
from scipy.optimize import curve_fit
from models import Models
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
    file_name = "mol_smiles.csv"
    df = pd.read_csv(file_name, delimiter=",", encoding='utf-8')
    parameters, covariance = curve_fit(Models.model_1, (df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), df['HSP_exp'], maxfev = 500000)

    fit_data = Models.model_2((df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), *parameters)
    SE = np.sqrt(np.diag(covariance))
    print(parameters, SE)
    