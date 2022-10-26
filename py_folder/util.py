import pandas as pd
import re
from collections import Counter
from zipfile import ZipFile
import os
import json

"""
CD-HIT-est clusters for concatenated resfinder and card db
"""


# get and prepare cluster file
def hash_map():
    cd_hit = open("C:/Users/diogo/Desktop/ARG_CR/Thesis/files/CD_HIT.sorted", "r")
    clust_n = []
    clust_x = []
    ref_seq = []
    lista_valores = []
    for line in cd_hit:
        if line[0] == ">":
            label = line.replace("\n", "")
            clust_x.append(label)
            lista_valores.append("$")
        if "*" in line:
            ref = re.findall(r'(>.*?) .', line)[0]
            ref_seq.append(ref)
        if line[0] != ">":
            ref = re.findall(r'(>.*?) .', line)[0]
            lista_valores.append(ref[0:-3])
    for c in range(len(clust_x)):
        clust_name = clust_x[c] + ref_seq[c]
        clust_n.append(clust_name[0:-3])
    indices = [index for index, element in enumerate(lista_valores) if element == "$"]
    lista_de_listas = []
    for l in range(len(indices) - 1):
        lista_de_listas.append(lista_valores[indices[l] + 1:indices[l + 1]])
    hash_map = dict(zip(clust_n, lista_de_listas))
    # fixing bug: cluster >Cluster 1012     0	903nt, > erm(N)_2_MZ015744... *
    hash_map[">Cluster 1012> erm(N)_2_MZ015744"] = hash_map.pop(">Cluster 10")
    hash_map[">Cluster 1012> erm(N)_2_MZ015744"] = ["> erm(N)_2_MZ015744"]
    return hash_map


# convert dic_c arg to cd-hit clusters
def clust_conversion(dic_c):
    converted = {}
    for id, l_arg in dic_c.items():
        for clust, l_seqs in hash_map().items():
            for arg in l_arg:
                if arg in l_seqs:
                    if id in converted.keys():
                        converted[id].append(clust)
                    else:
                        converted[id] = [clust]
    return converted


# read database
def read_db(database):
    with open(database, "r") as f:
        file = f.readlines()
        filer = [i.strip("\n") for i in file]
        lista_seqs = []
        seqs = ""
        for line in filer:
            if not line.startswith(">"):
                seqs += line
            else:
                lista_seqs.append(seqs + "\n")
                seqs = line + "\n"
    return lista_seqs


# writes fasta files from the clusters
def write_fasta(hash, lista_seqs_db, output_path):
    for file_names, ids_list in hash.items():
        filename = re.findall(">.*>", file_names)
        filename = re.sub(">", "", filename[0])
        filename = re.sub(" ", "", filename)
        with open(output_path + filename + ".fasta", "w") as wf:
            for ids in ids_list:
                for db_seq in lista_seqs_db:
                    if ids in db_seq:
                        wf.write(db_seq)
        wf.close()


# hash = hash_map()
# comb_db = read_db(database="C:/Users/diogo/Desktop/ARG_CR/Thesis/db/rf_card_db.fasta")
# write_fasta(hash, comb_db, "C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/clust_fasta/")


"""
Core Resistome
"""

# Create core with any prevalence percentage
def make_core(core_perc, comb_file, name):
    conta_i = 0
    conta_e = 0
    conta_s = 0
    conta_fw = 0
    conta_ww = 0
    for i in comb_file['type']:
        if i == "Influent":
            conta_i += 1
        if i == "Effluent":
            conta_e += 1
        if i == "Sludge":
            conta_s += 1
        if i == "Freshwater":
            conta_fw += 1
        if i == "Influent" or i == "Effluent":
            conta_ww += 1
    min_val_i = int(conta_i * core_perc)
    min_val_e = int(conta_e * core_perc)
    min_val_s = int(conta_s * core_perc)
    min_val_fw = int(conta_fw * core_perc)
    min_val_ww = int(conta_ww * core_perc)
    setted = []
    set_list = list(comb_file[name])
    for l in set_list:
        setted.append(list(dict.fromkeys(l)))
    type_arg = dict(zip(comb_file['id'], setted))
    id_type = dict(zip(comb_file['id'], comb_file['type']))
    estendido = {"Effluent": [], "Influent": [], "Sludge": [], "Freshwater": [], "Wastewater": []}
    for key, value in estendido.items():
        for id, arg in type_arg.items():
            if id_type[id] == key:
                for l in arg:
                    estendido[key].append(l)
        if key == "Effluent" or key == "Influent":
            for v in value:
                estendido["Wastewater"].append(v)
    final = estendido.copy()
    for key, value in estendido.items():
        if key == "Effluent":
            z = Counter(estendido[key])
            z = dict(z)
            myDict = {key: val for key, val in z.items() if val >= min_val_e}
            final[key] = myDict
        if key == "Influent":
            z = Counter(estendido[key])
            z = dict(z)
            myDict = {key: val for key, val in z.items() if val >= min_val_i}
            final[key] = myDict
        if key == "Sludge":
            z = Counter(estendido[key])
            z = dict(z)
            myDict = {key: val for key, val in z.items() if val >= min_val_s}
            final[key] = myDict
        if key == "Freshwater":
            z = Counter(estendido[key])
            z = dict(z)
            myDict = {key: val for key, val in z.items() if val >= min_val_fw}
            final[key] = myDict
        if key == "Wastewater":
            z = Counter(estendido[key])
            z = dict(z)
            myDict = {key: val for key, val in z.items() if val >= min_val_ww}
            final[key] = myDict
    return final


# Difference between soft core (75%) and full core (90%)
def soft_full(df_90, df_75):
    keep_me = {}
    temp = pd.DataFrame(index=df_90.index, columns=df_90.columns)
    soft = pd.DataFrame(index=df_90.index)
    for col in df_75:
        if col not in temp:
            soft[col] = df_75[col]
        else:
            for rowIndex, row in df_75.iterrows():
                for columnIndex, value in row.items():
                    for rI, r in df_90.iterrows():
                        for cI, v in r.items():
                            if columnIndex == cI and rowIndex == rI and value != 0 and v == 0:
                                keep_me[columnIndex] = {rowIndex: value}
    add = pd.DataFrame.from_dict(keep_me)
    soft = pd.concat([add, soft], axis=1)
    return soft


"""
Reads Relative Abundance
"""


# get ResFinder AMR class
def get_amr_class():
    amr_class = pd.read_csv('C:/Users/diogo/Desktop/ARG_CR/Thesis/files/phenotypes.txt', sep='\t')
    for line in amr_class["Gene_accession no."]:
        if "-" in line:
            info = re.findall(r'(.*?)-', line)
        else:
            info = re.findall(r'(.*?)_', line)
        amr_class = amr_class.replace(line, info[0])
    amr_dict = dict(zip(amr_class["Gene_accession no."], amr_class["Class"]))
    return amr_dict


# dic = {sra_id : n_reads}
def map_reads():
    directory = "C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/reads_count_unzip"
    n_reads = {}
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            g = open(f)
            data = json.load(g)
            for key, value in data.items():
                if key == "metadata":
                    for iter in value:
                        for x, y in iter[10].items():
                            if x == "read_count":
                                n_reads[iter[1]] = y
    return n_reads


"""
Cluster consensus
"""


def get_clust_hit(DIR, folder):
    id_clustlist = {}
    for id in lista_ids():
        caminho = DIR + folder + id + "/" + id + "_kma.res"
        con_df = pd.read_csv(caminho, sep='\t')
        con_df.rename({'#Template': id}, axis=1, inplace=True)
        id_clustlist[id] = list(con_df[id])
    return id_clustlist


def get_con_df(df, dicio):
    marked_df = df.copy()
    import numpy as np
    extended = []
    for key, value in dicio.items():
        for v in value:
            if v not in extended:
                extended.append(v)
            else:
                continue
    for i in extended:
        marked_df[i] = np.NaN
    for key, value in dicio.items():
        for v in value:
            if v in marked_df.columns and key in marked_df.index:
                marked_df.at[key, v] = 1
    return marked_df


"""
Taxonomy
"""


def unzip_files(dir_name, extraction):
    import os, zipfile
    extension = ".zip"
    os.chdir(dir_name)  # change directory from working dir to dir with files
    for item in os.listdir(dir_name):  # loop through items in dir
        if item.endswith(extension):  # check for ".zip" extension
            file_name = os.path.abspath(item)  # get full path of files
            zip_ref = zipfile.ZipFile(file_name)  # create zipfile object
            zip_ref.extractall(extraction)  # extract file to dir
            zip_ref.close()  # close file
            # os.remove(file_name)  # delete zipped file


"""
Other Utils
"""


# get sra id list
def lista_ids():
    file = "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/id_list.txt"
    handle = open(file, "r")
    lista_ids = handle.read().splitlines()
    return lista_ids

