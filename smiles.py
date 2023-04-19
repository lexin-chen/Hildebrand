import pandas as pd
from urllib.request import urlopen

def CIRconvert(ids):
    try:
        url = f'http://cactus.nci.nih.gov/chemical/structure/{ids}/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        print(f"{ids} not work")
        return 'Did not work'

if __name__ == '__main__':
    file_name = "all_molecules.csv"
    df = pd.read_csv(file_name, delimiter=",", encoding='utf-8')
    df['name'] = df['name'].astype(str)
    df['SMILES'] = df.iloc[:, 0].apply(CIRconvert)
    df.to_csv('mol_smiles.csv', index=False)
    