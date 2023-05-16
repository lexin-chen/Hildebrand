import os
import pandas as pd
import modules as mod

if __name__ == "__main__":
    """Plot the models for one type of dipoles 
    (all_dipoles, low_dipoles, high_dipoles)
    """
    dipole = "all_dipoles"
    save_dir_path = f"plots/{dipole}"
    if not os.path.exists(save_dir_path):
        os.makedirs(save_dir_path)
    text = " ".join(word.capitalize() for word in dipole.split("_"))
    df = pd.read_csv(f"{dipole}.csv", delimiter=",", encoding='utf-8')
    exp = df['HSP_exp']
        
    mod.loop_models_plot(df, exp, text, save_dir_path)

    """Plot the models for all types of dipoles
    """
    for dipole in ["all_dipoles", "low_dipoles", "high_dipoles"]:
        save_dir_path = f"plots/{dipole}"
        if not os.path.exists(save_dir_path):
            os.makedirs(save_dir_path)
        text = " ".join(word.capitalize() for word in dipole.split("_"))
        df = pd.read_csv(f"{dipole}.csv", delimiter=",", encoding='utf-8')
        exp = df['HSP_exp']
          
        mod.loop_models_plot(df, exp, text, save_dir_path)