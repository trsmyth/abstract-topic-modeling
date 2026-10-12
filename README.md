# abstract-topic-modeling
Abstract Topic Modeling: An NLP Alternative to Conventional Pathway Enrichment

Pathway enrichment is a common bioinformatic technique which attempts to identify biologically meaningful pathways demonstrating over-representation related to a list of genes or proteins. This enrichment generally follows differential expression analysis and relies on curation of pathway-gene/protein links. While these techniques, such as gene ontology (GO) enrichment, are mainstays of -omics analysis, they rely on public, expert input to create and update pathways. Additionally, in my experience, a large number of irrelevant pathways demonstrate enrichment, requiring significant amounts of time to filter out the massive noise generated during traditional analysis. As such, these techniques often restrict up- and downstream analysis, forcing a 1:1 DE approach followed by time-intensive pathway curation which is already restricted to a previously characterized, possibly outdated context.

This repository contains code describing a topic modeling alternative to conventional pathway enrichment. Following the desired form of gene/protein enrichment, such as a differential expression, clustering, or regression, Pubmed title/abstract text data is accessed through the Entrez API. This text data is linked to their related genes/proteins and assessed for topic enrichment. This approach allows for identification of emerging trends or the identification of biologically meaningful topics which were previously split due to investigator discoveries between cell types, tissues, or fields.

- project_functions.py: Functions related to the project. Includes functions for the linking of current gene names to previous symbols, querying Entrez for text data download, and assembly of downloaded text data into a single_gene class.
- Gene_db.py: Script for the creation of a local SQL database for the storage of Pubmed text data.
- Entrez_query.py: Script for Pubmed query and text download. Script allows for tracking of previously downloaded PMIDs to prevent redundant queries and downloads.
- Topic_modeling.py: Script for initial approach toward BERTopic topic modeling. This includes loading and sorting of text data, topic modeling, collapse and combination of related keywords, and clustering/visualization of identified topics.

Ongoing development is related to topic modeling options. The current approach utilizes BERTopic to identify distinct clusters of related text and corresponding key terms. Related terms are combined and overarching topics are identified using clustering of text embedding.

Below is an early example. Genes enriched in airway basal cells were determined using elastic net and clustered based on Spearman correlation. Genes within cluster 4 were ranked by elastic net correlations and gene set enrichment analysis (GSEA) of GO terms were assessed. 

<img height="150" alt="image" src="https://github.com/user-attachments/assets/1a8d65f2-a799-4ce5-853e-9a4375673f19" />

Overall, very little enrichment is detected. Additionally, terms are incredibly broad, suggesting no unifying pathway enrichment. In contrast, below is a plot demonstrating initial topic enrichment.

<p align="center">
  <img height="400" alt="image" src="https://github.com/user-attachments/assets/f705c1d0-ce13-424a-812a-988e8cc2c8f9" />
</p>

While smaller gene sets are expected to enrich generic terms such as those found in cluster 0, enrichment of terms related to oxidase enzymes, blood, and the liver introduce a potential starting point for further analysis.

This project is ongoing and will attempt to incorporate further gene, pathway, and topic enrichment techniques.
