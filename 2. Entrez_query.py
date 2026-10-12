import time
import sqlite3
import pandas as pd
from project_functions import query_entrez, create_seed_terms

directory = './Input'

# Example dataset from previous run of an ongoing project.
# This data will be uploaded when the project related project is submitted for publication.
Cluster_genes = pd.read_csv(f'{directory}/All_Cluster_Results.csv')

gene_names = Cluster_genes.loc[Cluster_genes.Max_coef_names == 'Basal', 'Gene_symbol'].tolist()

#######################################################

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

#######################################################

api_key = ""
email = ""

query = """SELECT PMID FROM pubmed_paper;"""
cursor.execute(query)
papers_in_db = [str(row[0]) for row in cursor.fetchall()]

query = """SELECT gene_name FROM gene_pmid_map;"""
cursor.execute(query)
genes_in_db = [str(row[0]) for row in cursor.fetchall()]

# Exists for initial db build. Can be replaced with: 
# genes_to_query = list(set(gene_names))
# if there is a desire to check for new manuscripts
genes_to_query = list(set(gene_names) - set(genes_in_db))
print(f'{len(genes_to_query)} genes left to query.')

#######################################################

for i in range(0, len(genes_to_query), 10):

    batch_genes = genes_to_query[i : (i + 10)]

    for gene in batch_genes:

        entrez_gene_data, new_pmids, all_pmids = query_entrez(gene, 
                                                              gene_term_dict = seed_data,
                                                              api_key = api_key, 
                                                              email = email, 
                                                              pmid_track = papers_in_db)
                
        papers_in_db.extend(new_pmids)

        ###########

        insert_query = '''
        INSERT OR IGNORE INTO pubmed_paper (pmid, title, abstract, keywords, mesh)
        VALUES (?, ?, ?, ?, ?);
        '''
        pmid_data = [tuple(d.values()) for d in entrez_gene_data]
        cursor.executemany(insert_query, pmid_data)

        ###########

        insert_query = '''
        INSERT OR IGNORE INTO gene_pmid_map (gene_name, pmid) 
        VALUES (?, ?);
        '''
        result = tuple((gene, item) for item in all_pmids)
        cursor.executemany(insert_query, result)
        conn.commit()

        print(f'{len(new_pmids)} of {len(all_pmids)} PMIDs added for gene {gene}.')

    print(f'Completed genes {i + 10} of {len(genes_to_query)}.')
    time.sleep(5)
