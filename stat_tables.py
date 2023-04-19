import os
import pandas as pd
import modules as mod

if __name__ == "__main__":
    save_dir_path = f"tables"
    for dipole in ["all_dipoles", "low_dipoles", "high_dipoles"]:
        if not os.path.exists(save_dir_path):
            os.makedirs(save_dir_path)
        text = " ".join(word.capitalize() for word in dipole.split("_"))
        df = pd.read_csv(f"{dipole}.csv", delimiter=",", encoding='utf-8')
        exp = df['HSP_exp']
          
        mod.loop_model_table(df, exp, save_dir_path, dipole, text)