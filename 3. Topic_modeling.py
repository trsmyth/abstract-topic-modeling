from sentence_transformers import SentenceTransformer
from sentence_transformers.models import StaticEmbedding

static = StaticEmbedding.from_model2vec("neuml/pubmedbert-base-embeddings-8M")
biobert_model = SentenceTransformer(modules = [static])

########################################################################

import sqlite3
import numpy as np
import pandas as pd
from project_functions import Single_gene, create_seed_terms, determine_distance, get_shingles

conn = sqlite3.connect('gene_cluster_abstract_data.db')
cursor = conn.cursor()

cursor.execute("PRAGMA table_info(gene_information)")
hugo_columns = [row[1] for row in cursor.fetchall()]

query = """SELECT * FROM gene_information"""
cursor.execute(query)
results = cursor.fetchall()

HUGO_dataset = pd.DataFrame.from_records(results, columns = hugo_columns)
Columns_of_interest = ['symbol', 'gene_name', 'alias_symbol', 'alias_name', 'prev_symbol', 'prev_name', 'ensembl_gene_id']

seed_data, seed_data_values = create_seed_terms(HUGO_dataset.loc[:, Columns_of_interest])

Columns_of_interest = ['symbol', 'alias_symbol', 'prev_symbol', 'ensembl_gene_id']
limited_seed_data, limited_seed_data_values = create_seed_terms(HUGO_dataset.loc[:, Columns_of_interest])

########################################################################

Cluster_genes = pd.read_csv(f'./Input/All_Cluster_Results.csv')

basal = Cluster_genes.loc[Cluster_genes.Max_coef_names == 'Basal', ]

cluster = 4
cluster_genes = basal.loc[basal.Cluster == cluster, 'Gene_symbol'].tolist()
gene_descriptions = basal.loc[basal.Cluster == cluster, 'Gene_summary'].tolist()

primary_query = """
SELECT pubmed_paper.pmid, pubmed_paper.title, pubmed_paper.abstract, pubmed_paper.keywords, pubmed_paper.mesh
FROM gene_pmid_map 
JOIN pubmed_paper USING (pmid)
WHERE gene_pmid_map.gene_name = ?;"""

all_text = []
text_number = {}
all_gene_seed_data = []

for i in range(len(cluster_genes)):
    gene = cluster_genes[i]
    description = gene_descriptions[i]

    cursor.execute(primary_query, (gene, ))
    gene_pmid_data = cursor.fetchall()
    pmid_data = [dict(zip(['PMID', 'Title', 'Abstract'], item)) for item in gene_pmid_data]

    single_gene = Single_gene(gene, description)

    try:
        gene_seed_data = seed_data[gene]

    except:
        try:
            gene = [k for k, v in seed_data.items() if Single_gene.preprocess_text([gene])[0] in v][0]
            gene_seed_data = seed_data[gene]
        except:
            gene_seed_data = None

    if gene_seed_data != None:
        all_gene_seed_data.append(gene_seed_data)

    for j in range(len(pmid_data)):
        single_gene.add_PMID_data(pmid_data[j], search_terms = gene_seed_data)

    gene_text = single_gene.isolate_gene_text(text_type = ['Title', 'Abstract'])
    
    all_text.extend(gene_text)
    text_number[gene] = len(gene_text)

encoded_txt = biobert_model.encode(all_text, show_progress_bar = True)
text_encodings = np.vstack(encoded_txt)

targets = []
gene_target = []
for i, key in enumerate(text_number.keys()):
    values = text_number[key]
    targets.extend(np.repeat(i, values))
    gene_target.extend(np.repeat(key, values))

########################################################################

import math
from bertopic import BERTopic
from cuml.cluster import KMeans
# from cuml.cluster import HDBSCAN
from cuml.manifold.umap import UMAP
from bertopic.representation import KeyBERTInspired
from sklearn.feature_extraction.text import CountVectorizer

representation_model = KeyBERTInspired(random_state = 12345)
cluster_model = KMeans(n_clusters = math.ceil(math.log(len(cluster_genes))), random_state = 12345, init = 'k-means++')
# cluster_model = HDBSCAN(min_cluster_size = 15, metric = 'euclidean', cluster_selection_method = 'eom', prediction_data = True)
vectorizer_model = CountVectorizer(stop_words = "english", ngram_range = (1, 2))

umap_model = UMAP(n_neighbors = 10, n_components = 15,
                min_dist = 0.0, metric = "cosine",
                init = 'random', random_state = 12345)

topic_model = BERTopic(calculate_probabilities = False, 
                    embedding_model = biobert_model, 
                    hdbscan_model = cluster_model,
                    umap_model = umap_model, 
                    representation_model = representation_model,
                    seed_topic_list = all_gene_seed_data,
                    vectorizer_model = vectorizer_model)

topics, _ = topic_model.fit_transform(all_text, embeddings = text_encodings)

top_topics = topic_model.get_topic_info()
topic_words = [item for sublist in top_topics['Representation'] for item in sublist]
filtered_topics = [item for item in topic_words if 'gene' not in item.split()]
filtered_list = [item for item in filtered_topics if 'protein' not in item.split()]
filtered_list = [item for item in filtered_topics if all(word.upper() not in cluster_genes for word in item.split())]

extracted_genes = [item for item in filtered_list 
                    if any(word in limited_seed_data_values for word in item.split()) 
                    or item.replace(" ", "") in limited_seed_data_values]

extracted_topics = [item for item in filtered_list if item not in extracted_genes]

flat_topics = sorted(extracted_topics, key = len, reverse = True)
unique_flat_topics = np.unique(flat_topics)

####################

for i, word in enumerate(unique_flat_topics):
    for j, word2 in enumerate(unique_flat_topics):
        if j <= i:
            continue
        if word.replace(" ", "") == word2.replace(" ", ""):
            if len(word) > len(word2):
                flat_topics = [item.replace(word2, word) for item in flat_topics]
            else:
                flat_topics = [item.replace(word, word2) for item in flat_topics]

topic_tracker = {}
for i, word1 in enumerate(flat_topics):
    if word1 in [item for sublist in list(topic_tracker.values()) for item in sublist]:
        continue      

    if word1 not in topic_tracker.keys():
        topic_tracker[word1] = [word1]
    
    word_in_values = []
    for j, word2 in enumerate(flat_topics):
        if j <= i:
            continue

        if len(set(word1.split()).intersection(word2.split())) > 0:
            topic_tracker[word1].append(word2)
            continue

        if len(word1) < 3 or len(word2) < 3:
            continue

        shingled_word1 = get_shingles(word1, n = 3)
        shingled_word2 = get_shingles(word2, n = 3)
        intersection = len(shingled_word1.intersection(shingled_word2))
        union = len(shingled_word1.union(shingled_word2))
        sim_score = intersection / len(shingled_word1)
        if sim_score >= 0.5:
            topic_tracker[word1].append(word2)
            continue

        for word in topic_tracker[word1]:
            if len(set(word.split()).intersection(word2.split())) > 0:
                word_in_values.append(word2)
    
    topic_tracker[word1].extend(word_in_values)

all_items = [item for item in list(topic_tracker.values())]
all_topics = [[item + "<br/>" for item in sublist] for sublist in all_items]
all_topics = [''.join(item) for item in all_topics]

################################

from sklearn.preprocessing import normalize
from sklearn.cluster import AgglomerativeClustering

embeded_topics = np.vstack([np.array(biobert_model.encode(item)) for item in topic_tracker.keys()])
normalized_data = normalize(embeded_topics, norm = 'l2', axis = 1) 

optimal_k = math.ceil(math.sqrt(len(all_topics)))

final_clustering = AgglomerativeClustering(n_clusters = optimal_k, 
                                            metric = 'cosine', 
                                            linkage = 'complete')

final_labels = final_clustering.fit_predict(normalized_data)

cluster_results = {}
for i in np.unique(final_labels):
    index = np.where(final_labels == i)[0]
    clustered_topics = np.array(list(topic_tracker.keys()))[index]

    merged_topics = []
    for single_topic in clustered_topics:
        linked_topics = topic_tracker[single_topic]
        merged_topics.append(single_topic + ':' + '\n' + "\n".join(linked_topics))

    sorted_list = sorted(merged_topics, key = lambda x: -x.count('\n'))
    cluster_results[i] = sorted_list

df = pd.DataFrame({k: pd.Series(v) for k, v in cluster_results.items()})
df = df.fillna("")
df.to_csv(f'Cluster_{cluster}.csv')

################################

import mpld3
import pacmap
import alphashape
from mpld3 import plugins
import matplotlib.pyplot as plt
import matplotlib.patches as patches

cmap = ['blue', 'orange', 'green', 'red', 'purple', 
        'gray', 'yellowgreen', 'pink', 'black', 'slateblue', 
        'cornflowerblue', 'darksalmon', 'forestgreen', 'cyan', 'chocolate', 
        'gold', 'crimson', 'mediumpurple', 'teal', 'sandybrown', 
        'firebrick', 'darkkhaki', 'mediumseagreen', 'lightcoral', 'powderblue']

model = pacmap.PaCMAP(n_components = 2, 
                        n_neighbors = 5, 
                        MN_ratio = 0.5, 
                        FP_ratio = 2.0) 

topic_embeds = model.fit_transform(embeded_topics)

fig, ax = plt.subplots()
for i in np.unique(final_labels):
    index = np.where(final_labels == i)[0]
    data_subset = topic_embeds[index][:, 0:2]
    center_point = [data_subset[:, 0].mean(), data_subset[:, 1].mean()]

    distance = []
    for point in data_subset:
        distance.append(determine_distance(center_point, point))

    sort_distance = np.argsort(distance)
    ordered_data_subset = data_subset[sort_distance][0:math.ceil(len(distance)*0.95)] 
    ordered_data_subset = ordered_data_subset + np.random.normal(0, 1e-9, ordered_data_subset.shape)
    hull = alphashape.alphashape(ordered_data_subset, alpha = 0.95)
    
    color = cmap[i]
    plot = plt.scatter(data_subset[:, 0], data_subset[:, 1], c = color, label = f'Cluster {i}')
    tooltip = mpld3.plugins.PointHTMLTooltip(plot, labels = [all_topics[k] for k in index])
    plugins.connect(fig, tooltip)

    if hull.geom_type == 'GeometryCollection':
        continue

    if hull.geom_type == 'Point':
        continue

    if hull.geom_type == 'LineString':
        coords = list(hull.coords.xy)
        ax.plot(list(coords[0]), list(coords[1]), color = color, linewidth = 2)

    else:
        try:
            hull_pts = hull.exterior.coords.xy
            coords = np.column_stack((hull_pts[0].tolist(), hull_pts[1].tolist()))
            ax.add_patch(patches.Polygon(coords, facecolor = color, edgecolor = 'black', linewidth = 1, alpha = 0.5))
        except:
            geoms = list(hull.geoms)
            for geom in geoms:
                hull_pts = geom.exterior.coords.xy
                coords = np.column_stack((hull_pts[0].tolist(), hull_pts[1].tolist()))
                ax.add_patch(patches.Polygon(coords, facecolor = color, edgecolor = 'black', linewidth = 1, alpha = 0.5))

    plt.subplots_adjust(right = 0.75)
    plt.legend(bbox_to_anchor = (1.05, 0.5), loc = "center left")
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')

    html_str = mpld3.fig_to_html(fig)

    file = f"Cluster_{cluster}.html"
    with open(file, "w") as f:
        f.write(html_str)
