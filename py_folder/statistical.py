import sys
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
sys.path.append("C:/Users/diogo/Desktop/ARG_CR/Thesis/py_folder")
from util import *
from core_resistome import *
from collections import defaultdict
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


"""
Abundância Relativa das reads
"""

# get reads
dir_zip_2 = "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/read_count/"
dir_end = ".JSON.zip"
for id_2 in lista_ids():
    DIR_sample_zip_2 = dir_zip_2 + id_2 + dir_end
    try:
        with ZipFile(DIR_sample_zip_2, 'r') as zipObj:
            listOfFileNames_2 = zipObj.namelist()
            for fileName in listOfFileNames_2:
                if fileName.startswith('KBase'):
                    zipObj.extract(fileName, 'C:/Users/diogo/Desktop/ARG_CR/Thesis/temp/reads_count_unzip')
    except:
        continue

reads = map_reads()
comb_file['n_reads'] = comb_file['id'].map(reads)

amr_dict = get_amr_class()
# excel the file
comb_file.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/Resume_table.xlsx")

# resfinder and silva
DIR = "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/kma/"
all_id_rpkm = {}
all_id_ra = {}
arg_ratio_per_reads = {}
silva_reads_hit = {}
all_arg_reads = {}
for id in lista_ids():
    DIR_sample = DIR + "resfinder/" + id
    DIR_sample_s = DIR + "silva/" + id
    dfs = pd.read_csv(DIR_sample_s + "/" + id + "_kma.mapstat", skiprows=6, delimiter="\t")
    df = pd.read_csv(DIR_sample + "/" + id + "_kma.mapstat", skiprows=6, delimiter="\t")
    df = df.rename(columns={"# refSequence": id})
    dfs = dfs.rename(columns={"# refSequence": id})
    Total_silva = dfs['readCount'].sum()
    silva_reads_hit[id] = Total_silva
    for line in df[id]:
        if "-" in line:
            info = re.findall(r'(.*?)-', line)
        else:
            info = re.findall(r'(.*?)_', line)
        df = df.replace(line, info[0])
    # adicionar coluna RPKM (amr reads por reads totais)
    for key, value in reads.items():
        if key == id:
            n1 = float(value) / 1000000
            rpm_list = []
            for i in range(len(df.readCount)):
                rpm = df['readCount'][i] / n1
                rpm_list.append(rpm)
                rpkm_list = []
            for r in range(len(df['bpTotal'])):
                n2 = float(df['bpTotal'][r] / 1000)
                rpkm = rpm_list[r] / n2
                rpkm_list.append(rpkm)
    df['RPKM'] = rpkm_list
    # coluna abundancia relativa (arg reads por arg reads totais)
    total_arg = df['readCount'].sum()
    all_arg_reads[id] = total_arg
    ra_list = []
    for i in range(len(df.readCount)):
        ra = df['readCount'][i] / total_arg
        ra_list.append(ra)
    df['RA'] = ra_list
    # Abundancia relativa: total arg reads/ total reads (por amostra)
    for key, value in reads.items():
        if key == id:
            t_reads = float(value)
            arg_per_total_r = total_arg / t_reads
            arg_ratio_per_reads[key] = arg_per_total_r
    # criar novo df com RPKM merged de acordo com nome do gene
    arg_rpkm = list(zip(df[id], df.RPKM))
    d = defaultdict(float)
    for x, y in arg_rpkm:
        d[x] += float(y)
    all_id_rpkm[id] = d
    # criar novo df com RA merged de acordo com nome do gene
    arg_ra = list(zip(df[id], df.RA))
    d1 = defaultdict(float)
    for x, y in arg_ra:
        for key, value in amr_dict.items():
            if x == key:
                d1[value] += float(y)
    all_id_ra[id] = d1

merged_rpkm = pd.DataFrame(all_id_rpkm).T.fillna(0)
merged_ra = pd.DataFrame(all_id_ra).T.fillna(0)

# get top 10 and Other
id_type = dict(zip(comb_file['id'], comb_file['type']))
medias = merged_ra.median()
medias_dict = dict(medias)
selection_top10 = sorted(medias_dict, key=medias_dict.get, reverse=True)[:10]
top10 = merged_ra[np.intersect1d(merged_ra.columns, selection_top10)]
top10["Other"] = 1 - top10.sum(axis=1)
top10["type"] = top10.index.to_series().map(id_type)
first_column = top10.pop('type')
top10.insert(0, 'type', first_column)

# excel the top 10 amr
top10.to_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/Output/top10_amr_class.xlsx")

# bar plot relative abundance per sample
ax = merged_ra.plot.bar(stacked=True, width=0.4)
plt.title("Relative abundance of AMR Classes per Sample")
plt.xticks(rotation=90)
plt.show()

sns.set_context("paper", rc={"font.size": 25, "axes.titlesize": 16, "axes.labelsize": 14})
# Boxplot da Abundancia relativa: total arg reads/ total reads (por amostra)
comb_file['ARG reads (relative abundance)'] = comb_file['id'].map(arg_ratio_per_reads)
sns.boxplot(x='type', y='ARG reads (relative abundance)', data=comb_file, order=["Influent", "Sludge", "Effluent", "Freshwater"])\
    .set(title="ResFinder: ARG Reads per Total Sample Reads", xlabel="Water Type", ylim=(0, 0.03))
# comb_file.boxplot(column=['RA_ARG'], by='type')
plt.show()

sns.set_context("paper", rc={"font.size": 25, "axes.titlesize": 16, "axes.labelsize": 14})
# Assembly Silva
assembly_silva = pd.read_excel("C:/Users/diogo/Desktop/ARG_CR/Thesis/files/assembly_silva.xlsx")
assembly_silva['16S rRNA contigs / Total contigs'] = pd.to_numeric(assembly_silva['silva_75%']) / pd.to_numeric(assembly_silva['Contigs'])
sns.boxplot(x='Type', y='16S rRNA contigs / Total contigs', data=assembly_silva, order=["Influent", "Sludge", "Effluent", "FreshWater"])\
    .set(title="SILVA Contigs (75%) per Total Contigs", xlabel="Water Type", ylim=(0, 1))
# comb_file.boxplot(column=['RA_ARG'], by='type')
plt.show()

sns.set_context("paper", rc={"font.size": 25, "axes.titlesize": 16, "axes.labelsize": 14})
# Boxplot SILVA
#comb_file['silva_reads'] = comb_file['id'].map(silva_reads_hit)
#comb_file['n_reads'] = comb_file['id'].map(all_arg_reads)
comb_file['16S rRNA reads / Total reads'] = pd.to_numeric(comb_file['silva']) / pd.to_numeric(comb_file['n_reads'])
a = sns.boxplot(x='type', y='16S rRNA reads / Total reads', data=comb_file, order=["Influent", "Sludge", "Effluent", "Freshwater"])\
    .set(title="SILVA Reads per Total Reads", xlabel="Water Type")
# comb_file.boxplot(column=['RA_ARG'], by='type')
plt.show()

sns.set_context("paper", rc={"font.size": 25, "axes.titlesize": 16, "axes.labelsize": 14})
comb_file['silva_reads'] = comb_file['id'].map(silva_reads_hit)
comb_file['all_arg_reads'] = comb_file['id'].map(all_arg_reads)
comb_file['ARG reads / 16S rRNA reads'] = pd.to_numeric(comb_file['all_arg_reads']) / pd.to_numeric(comb_file['silva_reads'])
sns.boxplot(x='type', y='ARG reads / 16S rRNA reads', data=comb_file, order=["Influent", "Sludge", "Effluent", "Freshwater"])\
    .set(title="ResFinder: ARG Reads per SILVA Reads", xlabel="Water Type", ylim=(0, 10))
# comb_file.boxplot(column=['RA_ARG'], by='type')
plt.show()


"""
PCA & unsupervised K-means
"""
# resfinder matriz 01
tudodegenes = []
for row in range(len(comb_file['resfinder'])):
    rf_ls = comb_file['resfinder'][row]
    tudodegenes.append(rf_ls)
flat_tudodegenes = [x for xs in tudodegenes for x in xs]
final_tg = list(dict.fromkeys(flat_tudodegenes))

# add all zeros
for i in final_tg:
    comb_file.insert(loc=9, column=i, value=0)
# if hit = 1
final_resfinder = comb_file.set_index('id')
id_rf = dict(zip(comb_file.id, final_resfinder['resfinder']))
for key, value in id_rf.items():
    for i in value:
        for z in final_resfinder.columns[8:]:
            if i == z:
                final_resfinder.loc[key, z] += 1
            else:
                continue

# standardização dos dados
scaler = StandardScaler()
final_resfinder_scaled = scaler.fit_transform(final_resfinder[final_tg])

# pca
rf_PCA = PCA()
rf_PCA.fit(final_resfinder[final_tg])
rf_pca_reduced = rf_PCA.transform(final_resfinder[final_tg])
# print(rf_PCA.explained_variance_ratio_)
print("2 PCA:", sum(rf_PCA.explained_variance_ratio_[0:2]))

barg = plt.figure(figsize=(5, 5))
plt.bar(["PC1", "PC2"], rf_PCA.explained_variance_ratio_[0:2])
plt.xlabel("Principal Components")
plt.ylabel("Percentagem")
plt.title("Variância explicada pelos dois primeiros PC")
plt.show()

for x in final_resfinder["type"].unique():
    sp = comb_file.index[final_resfinder["type"] == x] - 1
    plt.plot(rf_pca_reduced[sp, 0], rf_pca_reduced[sp, 1], "o", label=x)
plt.title("ResFinder: ARGs PCA")
plt.xlabel("PC1: " + str(rf_PCA.explained_variance_ratio_[0])[0:4])
plt.ylabel("PC2: " + str(rf_PCA.explained_variance_ratio_[1])[0:4])
plt.legend(loc="lower right")
plt.show()

# k means clustering resfinder
from sklearn.cluster import KMeans

data = final_resfinder[final_tg]
pca = PCA(2)
df_1 = pca.fit_transform(data)

kmeans = KMeans(n_clusters=4, max_iter=1000)
label = kmeans.fit_predict(df_1)
centroids = kmeans.cluster_centers_
u_labels = np.unique(label)
for i in u_labels:
    plt.scatter(df_1[label == i, 0], df_1[label == i, 1], label=i)
plt.scatter(centroids[:, 0], centroids[:, 1], s=80, color='k')
plt.title("Resfinder kmeans")
plt.legend()
plt.show()

# optimal kmeans resfinder
Sum_of_squared_distances = []
K = range(1, 15)
for k in K:
    km = KMeans(n_clusters=k)
    km = km.fit(df_1)
    Sum_of_squared_distances.append(km.inertia_)

plt.plot(K, Sum_of_squared_distances, 'bx-')
plt.xlabel('k')
plt.ylabel('Sum_of_squared_distances')
plt.title('Elbow Method For Optimal k: Resfinder')
plt.show()


#heatmap
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas import DataFrame
import seaborn as sns

df = pd.read_excel("C:/Users/diogo/Desktop/Output/consensus/freshwater_hm.xlsx", index_col=0)
sns.heatmap(df, cmap='BuPu', linewidths=0.5, vmax=1).set(title="Genes Presence Heatmap")
plt.show()
