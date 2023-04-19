import pandas as pd
import numpy as np
from enum import Enum
from scipy.optimize import curve_fit
import modules.models as models
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from adjustText import adjust_text
import os

def sep_by_dipole(file_name, df=None):
    if not df:
        df = pd.read_csv(file_name, delimiter=",", encoding='utf-8')
    low_dipole = df[df["Dipole"] < 5]
    high_dipole = df[df["Dipole"] >= 5]

    low_dipole.to_csv("low_dipole.csv", index=False)
    high_dipole.to_csv("high_dipole.csv", index=False)

def structure_annotation(smiles_list, exp_list, fit_data, ax):
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    
    for mol, x, y in zip(mol_list, exp_list, fit_data):
        drawer = rdMolDraw2D.MolDraw2DSVG(250, 250)
        opts = Draw.DrawingOptions()
        opts.bondLineWidth = 20
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        image = Draw.MolToImage(mol, size=(250, 250))
        x_offset = 0.3 + (250/1000) # calculate x offset based on image size
        ax.imshow(image, extent=[x+x_offset-0.3, x+x_offset+0.3, y-0.3, y+0.3])
    return ax

def loop_models_plot(df, exp, new, save_dir_path):
    for i in range(1, 12):
        func = getattr(models, f"model_{i}")
        parameters, covariance = curve_fit(func, (df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), 
                                        exp, maxfev = 500000)
        fit_data = func((df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), 
                        *parameters)
        SE = np.sqrt(np.diag(covariance))
        fig, ax = plt.subplots()
        # ax = structure_annotation(df['SMILES'], df['HSP_exp'], fit_data, ax)
        x = np.linspace(min(df['HSP_exp']), max(df['HSP_exp']))
        ax.scatter(df['HSP_exp'], fit_data, marker='o', zorder=10, color='#0047ab')
        ax.plot(x, x, dashes=[4, 4], color='#c57300')
        texts = []
        for j, txt in enumerate(df["formula"]):
            texts.append(ax.text(exp[j], fit_data[j], txt, ha='center', va='center'))
        adjust_text(texts, ax=ax)
        ax.title.set_text(f"Model {i}: {new}")
        ax.set_xlabel("Experimental Values")
        ax.set_ylabel("Theoretical Values")
        figure_name = f"model{i}"
        fig.savefig(f"{save_dir_path}/{figure_name}", bbox_inches = "tight", dpi = 500, transparent = True)
        print(f"Finished with Model {i}...")
    
def print_stats(formula, parameters, SE, mean_abs_err, rmse, loocv, aicc):
    """ Prints statistical results for a given formula, including mean absolute error (MAE), root-mean-square error (RMSE),
    leave-one-out cross-validation (LOOCV), and difference in Akaike information criterion corrected(dAICc).
    Args:
        formula (str): The formula used to fit the data.
        parameters (array-like): Array of the fitted parameters from `curve_fit`.
        SE (array-like): Array of the standard errors of the fitted parameters.
        mean_abs_err (float): The mean absolute error of the fitted data.
        rmse (float): The root-mean-square error of the fitted data.
        loocv (float): The leave-one-out cross-validation error of the fitted data.
        aicc (float): The Akaike information criterion corrected for finite sample sizes.
    Returns:
        stats_df (pandas.DataFrame): A DataFrame of the calculated statistics for the given formula.
    """   
    stats_dict = {
        "Model 1": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "γ"],
        "Model 2": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "ξ"],
        "Model 3": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "γ"],
        "Model 4": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "h"],
        "Model 5": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "a"],
        "Model 6": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "γ", "ξ", "h"],
        "Model 7": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "γ", "ξ", "a"],
        "Model 8": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b","ξ", "h"],
        "Model 9": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b" , "γ", "ξ", "a"],
        "Model 10": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "γ", "ξ", "h"],
        "Model 11": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "γ", "ξ", "a"],
    }   
    stats = {" ": stats_dict[formula],
            f"{formula}": [float('%.4g' % mean_abs_err), float('%.4g' % rmse), float('%.4g' % loocv), 
            float('%.4g' % aicc)] + [float('%.4g' % p) for p in parameters]}

    stats_df = pd.DataFrame(stats).set_index(" ")
    return stats_df

if __name__ == "__main__":
    dipole = "all_dipoles"
    new = " ".join(word.capitalize() for word in dipole.split("_"))
    save_dir_path = f"plots/{dipole}"
    if not os.path.exists(save_dir_path):
        os.makedirs(save_dir_path)
    df = pd.read_csv(f"{dipole}.csv", delimiter=",", encoding='utf-8')
    exp = df['HSP_exp']
    
    loop_models_plot(df, exp, new, save_dir_path)
    