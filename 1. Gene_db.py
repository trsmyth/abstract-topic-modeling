import sqlite3
import pandas as pd

# https://www.genenames.org/download/archive/
hugo_data = pd.read_csv('./Input/hgnc_complete_set.csv', dtype = str)

columns = ['hgnc_id', 'symbol', 'name', 'alias_symbol', 'alias_name', 'prev_symbol', 'prev_name', 
           'gene_group', 'gene_group_id', 'entrez_id', 'ensembl_gene_id', 'locus_group', 'location_sortable']

hugo_data = hugo_data.loc[:, columns].to_dict(orient = 'records')
hugo_data_db = [tuple(d.values()) for d in hugo_data]

########################################################

conn = sqlite3.connect('gene_cluster_abstract_data.db')
cursor = conn.cursor()
cursor.execute('PRAGMA foreign_keys = ON;')

cursor.execute("""
               CREATE TABLE IF NOT EXISTS gene_information (
               hgnc_id TEXT PRIMARY KEY,
               symbol TEXT,
               gene_name TEXT,
               alias_symbol TEXT,
               alias_name TEXT,
               prev_symbol TEXT,
               prev_name TEXT,
               gene_group TEXT,
               gene_group_id TEXT,
               entrez_id TEXT,
               ensembl_gene_id TEXT,
               locus_group TEXT,
               location_sortable TEXT
               )
               """)

conn.commit()

insert_query = '''
    INSERT OR IGNORE INTO gene_information (hgnc_id, symbol, gene_name, alias_symbol, alias_name, prev_symbol, 
    prev_name, gene_group, gene_group_id, entrez_id, ensembl_gene_id, locus_group, location_sortable) 
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    '''

cursor.executemany(insert_query, hugo_data_db)
conn.commit()

####################################################################

cursor.execute(
"""
CREATE TABLE IF NOT EXISTS pubmed_paper (
    pmid INTEGER PRIMARY KEY, 
    title TEXT, 
    abstract TEXT, 
    keywords TEXT, 
    mesh TEXT
);
"""
)

cursor.execute(
"""
CREATE TABLE IF NOT EXISTS gene_pmid_map (
    gene_name VARCHAR(50) NOT NULL,
    pmid INT REFERENCES pubmed_paper(pmid),
    PRIMARY KEY (gene_name, pmid)
);

"""
)

conn.commit()

'''
1. TABLE pubmed_paper is the manuscript data table. The manuscript PMID is the primary
key with manuscript text data saved in indicated columns. This will prevent the same
manuscript from being queried and saved multiple times.

2. TABLE gene_pmid_map acts as the intermediate between genes to be assessed for topics and pmid data.
The gene name will be a generic list of genes, for example, those linked to gene sets or DE results 
while pmid references the pmid in pubmed_paper.
'''

####################################################################

cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = [row[0] for row in cursor.fetchall()]

for table in tables:
    print(f"\nTable: {table}")
    cursor.execute(f"PRAGMA table_info({table});")
    
    for col in cursor.fetchall():
        print(f"--Column: {col[1]}")
