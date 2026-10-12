# abstract-topic-modeling
Abstract Topic Modeling: An NLP Alternative to Conventional Pathway Enrichment

Pathway enrichment is a common bioinformatic technique which attempts to identify biologically meaningful pathways demonstrating over-representation related to a list of genes or proteins. This enrichment generally follows differential expression analysis and relies on curation of pathway-gene/protein links. While these techniques, such as gene ontology (GO) enrichment, are mainstays of -omics analysis, they rely on public, expert input to create and update pathways. Additionally, in my experience, a large number of irrelevant pathways demonstrate enrichment, requiring significant amounts of time to filter out the massive noise generated during traditional analysis. As such, these techniques often restrict up- and downstream analysis, forcing a 1:1 DE approach followed by time-intensive pathway curation which is already restricted to a previously characterized, possibly outdated context.

This repository contains code describing a topic modeling alternative to conventional pathway enrichment. Following the desired form of gene/protein enrichment, such as a differential expression, clustering, or regression, Pubmed title/abstract text data is accessed through the Entrez API. This text data is linked to their related genes/proteins and assessed for topic enrichment. This approach allows for identification of emerging trends or the identification of biologically meaningful topics which were previously split due to investigator discoveries between cell types, tissues, or fields.

- project_functions.py: Functions related to the project. Includes functions for the linking of current gene names to previous symbols, querying Entrez for text data download, and assembly of downloaded text data into a single_gene class.
- Gene_db.py: Script for the creation of a local SQL database for the storage of Pubmed text data.
- Entrez_query.py: Script for Pubmed query and text download. Script allows for tracking of previously downloaded PMIDs to prevent redundant queries and downloads.

Ongoing development is related to topic modeling options. The current approach utilizes BERTopic to identify distinct clusters of related text and corresponding key terms. Related terms are combined and overarching topics are identified using clustering of text embedding.
