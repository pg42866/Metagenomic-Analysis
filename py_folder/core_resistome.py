"""
Import dos ficheiros, tratamento dos dados e construção dos df
"""
import sys

import numpy as np
import pandas as pd

sys.path.append("C:/Users/diogo/Desktop/ARG_CR/Thesis/py_folder")
from util import *

# read_template as df
comb_file = pd.read_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/files/final_dataset.xlsx")

"""
ALL ARG SET UP - KMA
"""
DIR = "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/kma/"

# resfinder
ficheiro_final = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_resfinder_list.txt", "w")
for id in lista_ids():
    DIR_sample = DIR + "resfinder/" + id
    try:
        ficheiro = open(DIR_sample + "/" + id + "_kma.fsa")
    except:
        continue
    ficheiro_final.write("$" + id + "\n")
    for line in ficheiro:
        ficheiro_final.write(line)
    ficheiro_final.write("\n")
    ficheiro.close()
ficheiro_final.close()

file = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_rf_list.txt")
dic_c = {}
for line in file:
    if line[0] == "$":
        dic_c[line[1:-1]] = []
        saver = line[1:-1]
    if line[0] == ">":
        dic_c[saver].append(line[:-1])

comb_file['resfinder'] = comb_file['id'].map(dic_c)
comb_file['resfinder_clust'] = comb_file['id'].map(clust_conversion(dic_c))

# card
ficheiro_final = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_card_list.txt", "w")
for id in lista_ids():
    DIR_sample = DIR + "card/" + id
    try:
        ficheiro = open(DIR_sample + "/" + id + "_kma.fsa")
    except:
        continue
    ficheiro_final.write("$" + id + "\n")
    for line in ficheiro:
        ficheiro_final.write(line)
    ficheiro_final.write("\n")
    ficheiro.close()
ficheiro_final.close()

file = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_card_list.txt")
dic_c = {}
for line in file:
    if line[0] == "$":
        dic_c[line[1:-1]] = []
        saver = line[1:-1]
    if line[0] == ">":
        dic_c[saver].append(re.findall(r'(.*?) \[', line)[0])

comb_file['card'] = comb_file['id'].map(dic_c)
comb_file['card_clust'] = comb_file['id'].map(clust_conversion(dic_c))

# silva hits
ficheiro_final = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_silva_list.txt", "w")
for id in lista_ids():
    DIR_sample = DIR + "silva/" + id
    try:
        ficheiro = open(DIR_sample + "/" + id + "_kma.fsa")
    except:
        continue
    ficheiro_final.write("$" + id + "\n")
    for line in ficheiro:
        ficheiro_final.write(line)
    ficheiro_final.write("\n")
    ficheiro.close()
ficheiro_final.close()

file = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/kma_silva_list.txt")
dic = {}
dicio = {}
for line in file:
    if line[0] == "$":
        dic[line[1:-1]] = []
        dicio[line[1:-1]] = None
        saver = line[1:-1]
    if line[0] == ">":
        info = [re.findall(r'(.*?) ', line)[0]]
        listToStr = ' '.join([str(elem) for elem in info])
        dic[saver].append(info)
for k, v in dic.items():
    dicio[k] = len(v)

comb_file['silva'] = comb_file['id'].map(dicio)

"""
Core Resistomes
"""
# Full core
# core_perc = 0.9
# cores
resfinder_core = make_core(0.9, comb_file, 'resfinder')
resfinder_clust_core = make_core(0.9, comb_file, 'resfinder_clust')
card_core = make_core(0.9, comb_file, 'card')
card_clust_core = make_core(0.9, comb_file, 'card_clust')
# df
df_resfinder_core = pd.DataFrame.from_dict(resfinder_core).fillna(0).T
df_card_core = pd.DataFrame.from_dict(card_core).fillna(0).T
df_resfinder_clust_core = pd.DataFrame.from_dict(resfinder_clust_core).fillna(0).T
df_card_clust_core = pd.DataFrame.from_dict(card_clust_core).fillna(0).T

# Soft core
# core_perc = 0.75
# cores
resfinder_soft_core = make_core(0.75, comb_file, 'resfinder')
resfinder_soft_clust_core = make_core(0.75, comb_file, 'resfinder_clust')
card_soft_core = make_core(0.75, comb_file, 'card')
card_soft_clust_core = make_core(0.75, comb_file, 'card_clust')
# df
df_soft_resfinder_clust_core = pd.DataFrame.from_dict(resfinder_soft_clust_core).fillna(0).T
df_soft_card_clust_core = pd.DataFrame.from_dict(card_soft_clust_core).fillna(0).T

# diferença entre soft e full
dif_core_resfinder = soft_full(df_resfinder_clust_core, df_soft_resfinder_clust_core).fillna(0)
dif_core_card = soft_full(df_card_clust_core, df_soft_card_clust_core).fillna(0)

# merged clust
rf_clust = dict(zip(comb_file.id, comb_file.resfinder_clust))
card_clust = dict(zip(comb_file.id, comb_file.card_clust))

# interseção
intersect = {}
for key, value in card_clust.items():
    for x, y in rf_clust.items():
        if key == x:
            if key not in intersect.keys():
                intersect[key] = []
            for i in value:
                for e in y:
                    if i == e:
                        n1 = value.count(i)
                        n2 = y.count(e)
                        minimo = min(n1, n2)
                        if i not in intersect[key]:
                            for z in range(minimo):
                                intersect[key].append(i)

comb_file['intersection_clust'] = comb_file['id'].map(intersect)
intersect_core = make_core(0.9, comb_file, 'intersection_clust')

# reuniao
reunion = {}
for key, value in card_clust.items():
    for x, y in rf_clust.items():
        if key == x:
            z = value + y
            if key not in reunion.keys():
                reunion[key] = []
                for i in z:
                    reunion[key].append(i)

comb_file['reunion_clust'] = comb_file['id'].map(reunion)
reunion_core = make_core(0.9, comb_file, 'reunion_clust')
df_reunion = pd.DataFrame.from_dict(reunion_core)

# excel the cores
df_resfinder_core.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/core/Core_rf.xlsx", sheet_name='resfinder_90%')
df_resfinder_clust_core.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/core/Core_rf_clust.xlsx",
                                 sheet_name='resfinder_clust_90%')
df_card_core.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/core/Core_card.xlsx", sheet_name='card_90%')
df_card_clust_core.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/core/Core_card_clust.xlsx",
                            sheet_name='card_clust_90%')

df_reunion.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/core/Core_reunion.xlsx")
cluster_number = []
for k, v in reunion_core.items():
    cluster_number.append(list(v.keys()))

cluster_number = [item for sublist in cluster_number for item in sublist]
cluster_number = list(dict.fromkeys(cluster_number))

"""
Consensus Seqs hits
"""

# hash = hash_map()
# comb_db = read_db(database="C:/Users/diogo/Desktop/ARG_CR/Thesis/db/rf_card_db.fasta")
# write_fasta(hash, comb_db, "C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/clust_fasta/")


con_file = pd.read_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/files/final_dataset.xlsx", index_col="id")
DIR = "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/kma_con/"

con_influent = get_clust_hit(DIR, "influent_con/")
con_sludge = get_clust_hit(DIR, "sludge_con/")
con_effluent = get_clust_hit(DIR, "effluent_con/")
con_freshwater = get_clust_hit(DIR, "freshwater_con/")
con_wastewater = get_clust_hit(DIR, "wastewater_con/")

df_con_influent = get_con_df(con_file, con_influent).fillna(0)
df_con_sludge = get_con_df(con_file, con_sludge).fillna(0)
df_con_effluent = get_con_df(con_file, con_effluent).fillna(0)
df_con_freshwater = get_con_df(con_file, con_freshwater).fillna(0)
df_con_wastewater = get_con_df(con_file, con_wastewater).fillna(0)

df_con_influent.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/consensus/df_con_influent.xlsx")
df_con_sludge.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/consensus/df_con_sludge.xlsx")
df_con_effluent.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/consensus/df_con_effluent.xlsx")
df_con_freshwater.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/consensus/df_con_freshwater.xlsx")
df_con_wastewater.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/consensus/df_con_wastewater.xlsx")

