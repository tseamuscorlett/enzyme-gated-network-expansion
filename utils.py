# packages
import csv
import ast
import copy
import pandas as pd
import matplotlib.pyplot as plt
import requests
from decimal import Decimal
from scipy.stats import spearmanr
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
from PIL import Image


# functions
def csv2dict(csv_file_path):
    result_dict = {}
    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        
        for row in csv_reader:
            try:
                result_dict[row[0]] = ast.literal_eval(row[1])  # value is list
            except:
                result_dict[row[0]] = row[1]
    return result_dict

def dict2csv(my_dict, csv_file_path):
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        for key, value in my_dict.items():
            csv_writer.writerow([key] + [value])

def todata(dict1, dict2, val_type = 'MEAN'):
    valid_keys = list(dict1.keys() & dict2.keys())
    data1 = [dict1[x] for x in valid_keys]
    data2 = [dict2[x] for x in valid_keys]
    
    if type(data1[0]) == dict:
        data1 = [x[val_type] for x in data1]
        
    if type(data2[0]) == dict:
        data2 = [x[val_type] for x in data2]
    
    return valid_keys, data1, data2

def spearman(dict1, dict2):
    valid_keys, data1, data2 = todata(dict1, dict2)
    correlation, p_value = spearmanr(data1, data2)
    p_value = '%.2E' % Decimal(p_value)
    return correlation, p_value

def scatter(dict1, dict2, x_axis = 'x-axis', y_axis = 'y-axis'):
    fig, ax = plt.subplots()
    
    valid_keys, data1, data2 = todata(dict1, dict2)
    plt.scatter(data1, data2, marker='o', color='b', alpha = 0.1, label='Data Points', zorder=2)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    # plt.savefig('scatter.png')
    ax.axline((0, 0), slope=1, color='k')
    plt.show()

def histogram(dict1, val_type = 'MEAN', bins = 10, x_axis = 'x-axis', y_axis ='counts'):
    data1 = list(dict1.values())
    
    if type(data1[0]) == dict:
        data1 = [x[val_type] for x in data1]
    
    plt.hist(data1, bins=bins, edgecolor='k')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.show()

def get_kegg_compound_smiles(kegg_id):
    url = f'http://rest.kegg.jp/get/{kegg_id}/mol'
    response = requests.get(url)
    if response.status_code == 200:
        mol_data = response.text
        mol = Chem.MolFromMolBlock(mol_data)
        return mol
    else:
        raise ValueError(f"Error fetching data for {kegg_id}: {cpd2name.get(kegg_id, '')}")

def draw_multiple_molecules(molecules, mol_labels=None):
    img = Draw.MolsToGridImage(molecules, molsPerRow=10, subImgSize=(400, 400), legends=mol_labels, maxMols=200)
    display(img)

    with open('cpds.png', mode='wb') as f:
        f.write(img.data)
        
def drawMols(molecule_kegg_ids):
    molecules = []
    labels = []

    # from SMILES
    for kegg_id in molecule_kegg_ids:
        try:
            mol_kegg = get_kegg_compound_smiles(kegg_id)
            if mol_kegg:
                molecules.append(mol_kegg)
                labels.append(cpd2name[kegg_id])
        except ValueError as e:
            print(e)
            
    if molecules:
        draw_multiple_molecules(molecules, labels)


# dictionaries
cpd2name = csv2dict('../data/assets/cpd2nameShort.csv')
cpd2nameLong = csv2dict('../data/assets/cpd2name.csv')

module2name = csv2dict('../data/assets/module2name.csv')
module2rns = csv2dict('../data/assets/module2rns.csv')

rn2cpds = csv2dict('../data/assets/rn2cpds_SI.csv')
rn2modules = csv2dict('../data/assets/rn2modules.csv')
rn2def = csv2dict('../data/assets/rn2def_versions.csv')
rn2rules = pd.read_pickle('../data/assets/rn2rules.20230224.pkl')
rn2eqn = csv2dict('../data/assets/rn2eqn_SI.csv')