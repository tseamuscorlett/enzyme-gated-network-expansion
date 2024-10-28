# packages
import csv
import ast
import copy
import requests
import pickle
import collections
import glob
import json
import random

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from decimal import Decimal
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import rankdata
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
from PIL import Image
from collections import Counter


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

def pearson(dict1, dict2):
    valid_keys, data1, data2 = todata(dict1, dict2)
    correlation, p_value = pearsonr(data1, data2)
    formatted_p_value = '{:e}'.format(p_value)
    return correlation, formatted_p_value

def scatter(dict1, dict2, x_axis = 'x-axis', y_axis = 'y-axis'):
    fig, ax = plt.subplots()
    
    valid_keys, data1, data2 = todata(dict1, dict2)
    plt.scatter(data1, data2, marker='o', color='b', alpha = 0.1, label='Data Points', zorder=2)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.show()

def histogram(dict1, val_type = 'MEAN', bins = 10, x_axis = 'x-axis', y_axis ='counts'):
    data1 = list(dict1.values())
    
    if type(data1[0]) == dict:
        data1 = [x[val_type] for x in data1]
    
    plt.hist(data1, bins=bins, edgecolor='k')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.show()

def loglog(dict1, dict2, x_axis = 'x-axis', y_axis = 'y-axis'):
    valid_keys, data1, data2 = todata(dict1, dict2)
    plt.scatter(np.log10(data1), np.log10(data2), marker='o', color='b', alpha = 0.1, label='Data Points', zorder=2)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.show()

def ylog(dict1, dict2, x_axis = 'x-axis', y_axis = 'y-axis'):
    valid_keys, data1, data2 = todata(dict1, dict2)
    plt.scatter(data1, np.log10(data2), marker='o', color='b', alpha = 0.1, label='Data Points', zorder=2)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    # plt.savefig('ylog_plot.svg', dpi=300, bbox_inches='tight')
    plt.show()
    
def xlog(dict1, dict2, x_axis = 'x-axis', y_axis = 'y-axis'):
    valid_keys, data1, data2 = todata(dict1, dict2)
    plt.scatter(np.log10(data1), data2, marker='o', color='b', alpha = 0.1, label='Data Points', zorder=2)
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

def rnWith(xgroup, rn2rules, rn2cpds):
    rnWithX = []
    for reaction, rules in rn2rules.items():
        if reaction in rn2cpds:  # only get rns that show up in SI (12872 -> 8558)
            for rule in rules:
                if xgroup in rule:
                    rnWithX.append(reaction)
                    break
    return rnWithX


# dictionaries
cpd2name = csv2dict('../data/assets/cpd2nameShort.csv')
cpd2nameLong = csv2dict('../data/assets/cpd2name.csv')

rn2cpds = csv2dict('../data/assets/rn2cpds_SI.csv')
rn2reac = csv2dict('../data/assets/rn2reac.csv')
rn2prod = csv2dict('../data/assets/rn2prod.csv')
rn2direction = csv2dict('../data/assets/rn2direction.csv')
rn2rev = csv2dict('../data/assets/rn2reversible.csv')
rn2modules = csv2dict('../data/assets/rn2modules.csv')
rn2def = csv2dict('../data/assets/rn2def_versions.csv')
rn2rules = pd.read_pickle('../data/assets/rn2rules.20230224.pkl')
rn2eqn = csv2dict('../data/assets/rn2eqn_SI.csv')

module2name = csv2dict('../data/assets/module2name.csv')
module2rns = csv2dict('../data/assets/module2rns.csv')

x2ns = csv2dict('../data/assets/xgroup2/xgroup2networkSize.csv')
x2rn = {}
for xgroup in x2ns.keys():
    x2rn[xgroup] = rnWith(xgroup, rn2rules, rn2cpds)
x2name = csv2dict('../data/assets/xgroup2/xgroup2name.csv')
x2class = csv2dict('../data/assets/xgroup2/xgroup2class.csv')
x2ds = csv2dict('../data/assets/xgroup2/DS_average/xgroup2DS_average_ArcBac_recovered.csv')
x2dsEuk = csv2dict('../data/assets/xgroup2/DS_average/xgroup2DS_eukaryotes_recovered.csv')
x2arc = csv2dict('../data/assets/xgroup2/xgroup2architecture.csv')
x2modules = csv2dict('../data/assets/xgroup2/xgroup2modules.csv')