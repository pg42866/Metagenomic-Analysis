"""
Taxonomy
"""
import sys
sys.path.append("C:/Users/diogo/Desktop/ARG_CR/Thesis/py_folder")
from core_resistome import *
from util import *
import glob

#unzip_files("C:/Users/diogo/Desktop/ARG_CR/Thesis/files/kaiju", "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/kaiju_unziped")

all_taxon = []
taxon_dict = {}
path = "C:/Users/diogo/Desktop/ARG_CR/Thesis/files/kaiju_unziped"
for filename in glob.glob(os.path.join(path, '*phylum.kaijuReport')):
    with open(os.path.join(os.getcwd(), filename), 'r') as f:
        dataf = pd.read_csv(filename, sep='\t')
        all_taxon.append(list(dataf['taxon_name']))
        for i in lista_ids():
            if i in filename:
                taxon_dict[i] = dict(zip(dataf['taxon_name'], dataf['percent']))

all_taxon = [item for sublist in all_taxon for item in sublist]
all_taxon = list(dict.fromkeys(all_taxon))
all_taxon.remove('unclassified')


taxon_df = pd.DataFrame.from_dict(taxon_dict).T.fillna(0)
taxon_df = taxon_df.drop(columns=['unclassified'])
taxon_df_binary = taxon_df.copy()
mapper = dict(zip(comb_file['id'], comb_file['type']))
taxon_df['type'] = taxon_df.index.map(mapper)
taxon_df_binary[taxon_df != 0] = 1
taxon_df_binary['type'] = taxon_df_binary.index.map(mapper)


# pca taxon_df
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 16})
taxon_PCA = PCA()
taxon_PCA.fit(taxon_df[all_taxon])
taxon_pca_reduced = taxon_PCA.transform(taxon_df[all_taxon])
print("2 PCA:", sum(taxon_PCA.explained_variance_ratio_[0:2]))

for x in taxon_df["type"].unique():
    sp = comb_file.index[taxon_df["type"] == x] - 1
    plt.plot(taxon_pca_reduced[sp, 0], taxon_pca_reduced[sp, 1], "o", label=x)
plt.title("Phylum Taxonomy PCA")
plt.xlabel("PC1: " + str(taxon_PCA.explained_variance_ratio_[0]))
plt.ylabel("PC2: " + str(taxon_PCA.explained_variance_ratio_[1]))
plt.legend(loc="lower right")
plt.show()


# pca taxon_df_binary
taxon_PCA = PCA()
taxon_PCA.fit(taxon_df_binary[all_taxon])
taxon_pca_reduced = taxon_PCA.transform(taxon_df_binary[all_taxon])
print("2 PCA:", sum(taxon_PCA.explained_variance_ratio_[0:2]))

for x in taxon_df_binary["type"].unique():
    sp = comb_file.index[taxon_df_binary["type"] == x] - 1
    plt.plot(taxon_pca_reduced[sp, 0], taxon_pca_reduced[sp, 1], "o", label=x)
plt.title("Phylum: Binary Taxonomy PCA")
plt.xlabel("PC1: " + str(taxon_PCA.explained_variance_ratio_[0]))
plt.ylabel("PC2: " + str(taxon_PCA.explained_variance_ratio_[1]))
plt.legend(loc="lower right")
plt.show()
