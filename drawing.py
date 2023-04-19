import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

# define the molecules
smiles_list = ['CCO', 'CNC', 'CCN']
mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
x_list = [1, 2, 3]
y_list = [1, 2, 3]

# create the plot
fig, ax = plt.subplots()

# plot each molecule
for mol, x, y in zip(mol_list, x_list, y_list):
    drawer = rdMolDraw2D.MolDraw2DSVG(250, 250)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    image = Draw.MolToImage(mol, size=(250, 250))
    ax.imshow(image, extent=[x-0.3, x+0.3, y-0.3, y+0.3])
plt.scatter(x_list, y_list, color='white', edgecolor='black', zorder=10)
# set the axis limits and labels
ax.set_xlim([0, 4])
ax.set_ylim([0, 4])
ax.set_xlabel('X')
ax.set_ylabel('Y')

# display the plot
plt.show()
