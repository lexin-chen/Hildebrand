import pandas as pd
import numpy as np
from enum import Enum
from scipy.optimize import curve_fit
import modules.models
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from adjustText import adjust_text
import os
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import LeaveOneOut, cross_val_score
from sklearn.linear_model import LinearRegression
from tabulate import tabulate

def sep_by_dipole(file_name, df=None):
    if not df:
        df = pd.read_csv(file_name, delimiter=",", encoding='utf-8')
    low_dipole = df[df["Dipole"] < 5]
    high_dipole = df[df["Dipole"] >= 5]

    low_dipole.to_csv("low_dipoles.csv", index=False)
    high_dipole.to_csv("high_dipoles.csv", index=False)

def loop_models_plot(df, exp, text, save_dir_path):
    for i in range(1, 12):
        func = getattr(modules.models, f"model_{i}")
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
        ax.title.set_text(f"Model {i}: {text}")
        ax.set_xlabel("Experimental Values")
        ax.set_ylabel("Theoretical Values")
        figure_name = f"model_{i}"
        fig.savefig(f"{save_dir_path}/{figure_name}", bbox_inches = "tight", dpi = 500, transparent = True)
        print(f"Finished with Model {i}...")

def loop_model_table(df, exp, save_dir_path, dipole, text):
    aicc_list = []
    stats_list = []
    for i in range(1, 12):
        func = getattr(modules.models, f"model_{i}")
        model = f"Model {i}"
        parameters, covariance = curve_fit(func, (df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), 
                                        exp, maxfev = 500000)
        fit_data = func((df['HOMO'], df['LUMO'], df['V_m'], df['chi'], df['eta']), 
                        *parameters)
        SE = np.sqrt(np.diag(covariance))
        mean_abs_err, rmse, loocv, aicc = get_errors(exp, fit_data, parameters)
        aicc_list.append(aicc)
        stats_df = print_stats(model, parameters, SE, mean_abs_err, rmse, loocv, aicc)
        stats_list.append(stats_df)
    stats_result = pd.concat(stats_list, axis=1)
    stats_result = replace_aicc(aicc_list, stats_result)
    stats_result = stats_result.fillna("-")
    with open(f"{save_dir_path}/{dipole}.txt", "w", encoding="utf-8") as f:
                    f.write(f"{text} Models\n")
                    f.write(tabulate(stats_result, headers="keys", stralign="center", tablefmt='fancy_grid'))

def get_errors(exp, fit_data, parameters):
    """ Calculates the mean absolute error (MAE), root-mean square error (RMSE), Akaike Information Criterion with correction for small samples (AICC),
    and RMSE of leave-one-out cross-validation (LOOCV) for a given property and its corresponding fit data.
    Args:
        exp (pandas Series): A 1-dimensional numpy array containing the actual property values.
        fit_data (np.ndarray): A 1-dimensional numpy array containing the fit data.
        parameters (np.ndarray): A 1-dimensional numpy array containing the parameters used for fitting from `curve_fit`.
    Returns:
        Tuple [float, float, float, float]: A tuple containing the MAE, RMSE, LOOCV RMSE and AICC, respectively.
    """
    mean_abs_err = mean_absolute_error(exp, fit_data)
    rmse = np.sqrt(mean_squared_error(exp, fit_data))
    aicc = calculate_aicc(exp.shape[0], rmse ** 2, parameters.shape[0]) 
    X = pd.DataFrame(fit_data)
    Y= pd.DataFrame(exp)
    cv = LeaveOneOut()
    model = LinearRegression()  
    scores = cross_val_score(model, X, Y, scoring='neg_mean_squared_error', cv=cv, n_jobs=-1)
    loocv = np.sqrt(np.mean(np.absolute(scores)))
    return mean_abs_err, rmse, loocv, aicc

def calculate_aicc(n, mse, num_params):
    """Calculates the Akaike information criterion (AIC) with small sample correction

    Args:
    n (int): number of observations
    mse (float): mean square error
    num_params (int): number of parameters

    Returns:
    aic_c (float): value for AIC corrected for small samples
    """
    aic = n * np.log(mse) + 2 * num_params
    aic_c = aic + (2 * num_params * (num_params + 1))/(n - num_params - 1)
    return (aic_c)

def replace_aicc(aicc_list, stats_result):
    """ Replace AICc values in stats_result DataFrame with the difference from the minimum AICc value.
    Args:
        aicc_list (list): AICc values for each regression model.
        stats_result (pandas.DataFrame): Results summary DataFrame with model statistics.
    Returns:
        stats_result (pandas.DataFrame): Updated DataFrame replaced with the AICc differences or dAICc.
    """
    if len(aicc_list) == 0:
        aicc_min = np.nan
    else:
        aicc_min = np.min(aicc_list)
    daicc_list = list(map(lambda v: v - aicc_min, aicc_list))
    stats_result.iloc[3]= ['%.2f' % elem for elem in daicc_list]
    return stats_result

def print_stats(model, parameters, SE, mean_abs_err, rmse, loocv, aicc):
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
        "Model 1": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b"],
        "Model 2": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b"],
        "Model 3": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b"],
        "Model 4": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "h"],
        "Model 5": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "a"],
        "Model 6": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "h"],
        "Model 7": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "a"],
        "Model 8": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "h"],
        "Model 9": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b" , "a"],
        "Model 10": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "h"],
        "Model 11": ["MAE", "RMSE", "LOOCV", "dAICc", "p", "b", "a"],
    }   
    stats = {" ": stats_dict[model],
            f"{model}": [float('%.4g' % mean_abs_err), float('%.4g' % rmse), float('%.4g' % loocv), 
            float('%.4g' % aicc)] + [float('%.4g' % p) for p in parameters]}

    stats_df = pd.DataFrame(stats).set_index(" ")
    return stats_df

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
        
if __name__ == "__main__":
    save_dir_path = f"tables"
    for dipole in ["all_dipoles", "low_dipoles", "high_dipoles"]:
        save_dir_path = f"plot/{dipole}"
        if not os.path.exists(save_dir_path):
            os.makedirs(save_dir_path)
        text = " ".join(word.capitalize() for word in dipole.split("_"))
        df = pd.read_csv(f"{dipole}.csv", delimiter=",", encoding='utf-8')
        exp = df['HSP_exp']
          
        #loop_models_plot(df, exp, text, save_dir_path)
        loop_model_table(df, exp)