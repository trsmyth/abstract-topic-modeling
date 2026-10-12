import math
import re
import json
import numpy as np
from Bio import Entrez
from langdetect import detect
from nltk.corpus import stopwords
import xml.etree.ElementTree as ET
from scipy.spatial.distance import cdist

########################################################################

def determine_distance(coordinate_1, coordinate_2):
    '''
    Given the coordinates of two points, determine the distance between them.
    '''    
    x0, y0 = coordinate_1
    x1, y1 = coordinate_2
    distance = math.sqrt((x1 - x0)**2 + (y1 - y0)**2)

    return distance

########################################################################

def get_shingles(text, n = 3):
    return {text[i:i+n] for i in range(len(text) - n + 1)}

########################################################################

def find_genes(unique_words, gene_dict, search_terms = None, strip_number_rear = 2, preprocessed = False):

    if not isinstance(gene_dict, dict):
        raise TypeError("gene_dict expects a dictionary with {'Gene' : [term]} structure.")
        
    if search_terms == None:
        search_terms = [values for subdict in gene_dict.values() for values in subdict]

    if preprocessed == True:
        search_terms = Single_gene.preprocess_text(search_terms)
        unique_words = Single_gene.preprocess_text(unique_words)

    overlapping_terms = []
    for word in unique_words:
        '''
        Check if the current word is in the search terms and save the match if present.
        
        If no match present, check if the word ends in a letter and number or number and letter. If yes, 
        from 1 to strip_number_rear + 1, strip the final i characters of woi and try again. This accounts for 
        situations with subunits which may have a '[0-9][a-z]' or '[a-z][0-9]' structure, i.e. XXXa1, XXXp1 
        or XXX1a which is not its own specified gene or gene product.
        '''
        overlapping_terms = []
        linked_keys = []

        for word in unique_words:
            if word in search_terms:
                matches = [term for term in search_terms if word == term]
                overlapping_terms.append((matches[0], word))
                linked_keys.append([k for k, v in gene_dict.items() if any(item in matches for item in v)])

            elif re.search(r'[a-z]\d+$', word) or re.search(r"\d[a-z]$", word):
                for i in range(1, strip_number_rear + 1): 
                    woi = word[:-i]
                    if len(woi) > 2 and (woi in search_terms):
                        matches = [term for term in search_terms if woi == term]
                        overlapping_terms.append((matches[0], word))
                        linked_keys.append([k for k, v in gene_dict.items() if any(item in matches for item in v)])
                        break

    return linked_keys, overlapping_terms

########################################################################

def query_entrez(gene, pmid_track = None, gene_term_dict = None, api_key = None, email = None, batch_size = 500):

    # Series of checks to ensure function parameters are correctly supplied
    # Prepares variables if optional parameters supplied
    if type(gene) != str:
        raise TypeError("gene expects a str containing the gene name.")
    
    if pmid_track != None and type(pmid_track) != list:
        raise TypeError("pmid_track expects a list of previously queried PMIDs.")
        
    if gene_term_dict != None and 'ENSG0' not in gene:
        if not isinstance(gene_term_dict, dict):
            raise TypeError("gene_term_dict expects a dictionary with {'Gene' : [term]} structure.")
        
        else:

            '''
            As HUGO symbols can and often do change over time, a previously published dataset could provide 
            an outdated symbol which has since been discontinued. As the gene_term_list is constructed from 
            HUGO provided datasets, this can cause future issues. Here, if no terms are found, it is probably 
            because the gene name is out of date of the gene key. In that case, check the dict values for the 
            gene name and use that list if found.
            '''
            try:
                values_list = gene_term_dict[gene]

            except:
                gene_key = [k for k, v in gene_term_dict.items() if Single_gene.preprocess_text([gene])[0] in v]
                values_list = gene_term_dict[gene_key[0]]
                if len(gene_key) > 1: # If the input gene is found in more than one primary HUGO ID, default to the first and warn the user.
                    print(f'Warning: Input gene {gene} is not a primary HUGO ID and matches {len(gene_key)} primary HUGO IDs. Defaulting to position 1: {gene_key[0]}')

            # Have ran into some problems with terms rather than words.
            # This puts every search term (and word) in quotes.
            search_terms = []
            for single_term in values_list:
                formatted_term = "\"" + single_term + "\"" + "[Abstract]"
                search_terms.append(formatted_term)

            # Join the search terms as a single query with OR between terms
            # Several examples of unrelated matches in titles or author blocks found in trials, 
            # so specify "[Abstract]" to specifically search the abstract text for the search terms
            query = " OR ".join(search_terms)

    else:
        query = gene

    if api_key != None:
        Entrez.api_key = api_key

    if email != None:
        Entrez.email = email

    ##########################################################################

    # Request data related to query and isolate the PMIDs
    handle = Entrez.esearch(db = "pubmed", term = query, retmax = 10000)
    results = Entrez.read(handle)
    handle.close()
    id_list = results['IdList']

    # If previously accessed/downloaded PMIDs are tracked, determine which 
    # PMIDs have not been accessed to access here
    if pmid_track != None:
        new_pmids = list(set(id_list) - set(pmid_track))

    else:
        new_pmids = batch_ids

    articles_data = [] # Storage list for extracted text data

    # Batching of requests is recommended by NCBI
    for i in range(0, len(new_pmids), batch_size):

        # Isolate PMIDs within current batch range
        batch_ids = new_pmids[i : (i + batch_size)]
        ids = ','.join(batch_ids)

        # Request PMID text data for current PMIDs
        handle = Entrez.efetch(db = "pubmed", retmode = "xml",  id = ids)
        xml_records = handle.read()
        handle.close()

        # Convert XML data to the root element of the parsed tree
        root = ET.fromstring(xml_records)

        # For each article, iterate over every element in the tree
        for article in root.iter('PubmedArticle'):

            # Try to find each element below
            try:

                # Extract PMID
                pmid = article.find('.//PMID').text

                #######################
                
                # Extract title
                title_text = ''
                title = article.find('.//ArticleTitle')

                # Some article titles have things like italics which breaks 
                # the normal text retrieval. This grabs the text and if
                # any sub elements are found, adds them and their tails.
                if title is not None:
                    title_text = title.text or ""
                    for sub in title:
                        if sub is not None:
                            title_text += sub.text or ""
                            title_text += sub.tail or ""

                #######################

                # Extract abstract
                abstract_text = ''

                # Abstracts frequently cause issues. This iterates over any
                # AbstractText elements and adds the found text, adding a space 
                # between each element which generally end after a period.
                for element in article.iter('AbstractText'):
                    if element is not None:
                        abstract_text += element.text or ""
                        for i, sub in enumerate(element):
                            if sub is not None:
                                abstract_text += sub.text or ""
                                abstract_text += sub.tail or ""
                                if i != len(element):
                                    abstract_text += ' '
                
                #######################

                # Extract keywords for article
                keywords_text = []

                # Again, keywords often have multiple subinstances which 
                # can cause problems. Here, iterate over the KeywordList 
                # and grab each keyword(s), adding each to a list.
                for element in article.iter('KeywordList'):
                    if element is not None:
                        for i, sub in enumerate(element):
                            if sub is not None:
                                keywords_text.append(sub.text)

                # Sometimes a None field is included in the keyterm list, possibly due to a 
                # formatting error when uploading/importing the keyword field to pubmed. This
                # removes None if found so it will play nice with join() below.
                keywords_text = [text for text in keywords_text if text != None]

                #######################

                # Extract MeSH terms
                mesh_terms = []

                for mesh_heading in article.iter('MeshHeading'):
                    descriptor = mesh_heading.find('DescriptorName').text
                    mesh_terms.append(descriptor)

                #######################

                text_data = {'PMID': pmid,
                             'Title': title_text,
                             'Abstract': abstract_text,
                             'Keywords' : ', '.join(keywords_text),
                             'MeSH Terms': ', '.join(mesh_terms)}
                
                articles_data.append(text_data)
            
            # If the element is not found, continue to the next one
            except AttributeError:
                continue

    # Return the article text data, list of new PMIDs extracted, and the total list of linked PMIDs
    return articles_data, new_pmids, id_list

########################################################################

def create_seed_terms(seed_dataframe):

    seed_data = {}

    for row_index in range(seed_dataframe.shape[0]):
        gene_seed_data = []
        row_gene_data_subset = seed_dataframe.iloc[row_index, :].dropna()
        gene_symbol = row_gene_data_subset.symbol
        row_gene_data_subset = row_gene_data_subset.drop('symbol')
        
        for column_index in range(len(row_gene_data_subset)):
            column_seed_data = []
            if "|" in row_gene_data_subset.iloc[column_index]:
                column_seed_data.extend(row_gene_data_subset.iloc[column_index].split('|'))
            else:
                column_seed_data.append(row_gene_data_subset.iloc[column_index])
            
            gene_seed_data.extend(column_seed_data)

        # Alternative gene names with only 3 characters have proven difficult, particularly 
        # if no numbers are present, as unrelated acronyms have frequently been confused 
        # for genes. Names with less than 4 characters are removed, and gene symbols are 
        # added after, as some official gene names are only 3 characters.
        
        gene_seed_data = [item for item in gene_seed_data if len(item) > 3]
        gene_seed_data.append(gene_symbol)
        gene_seed_data = Single_gene.preprocess_text(gene_seed_data)

        seed_data.update({gene_symbol : gene_seed_data})

    seed_data_values = [values for subdict in seed_data.values() for values in subdict]

    return seed_data, seed_data_values

########################################################################

def compute_medoid(points):
    
    dist_matrix = cdist(points, points, metric = 'euclidean') # Pairwise Euclidean distances
    medoid_idx = np.argmin(dist_matrix.sum(axis = 0)) # Index of the point with the minimum sum of distances
    return points[medoid_idx]

########################################################################

class Single_gene:
    '''
    Class containing data for a single gene. 
    Determines topics found within pubmed text data related to that gene which has been downloaded previously.

    Parameters:
    -----------
    gene_name: Name of the gene. Will be converted to lower case for ease of use later.

    description: Optional gene description. Generally will be NCBI description. Is added to other text data
    when assessing gene topics.

    bertopic_model: Optional BERTopic model (see bertopic package) which will be used for topic discovery. Can
    be added later on using 'set_bertopic' function.

    pmid_data: Names of fields for pmid data which will be added to the class. Defaults to PMID, Title, Abstract,
    Keywords, and MeSH Terms which are commonly found in PMID datasets. Can be modified to include whatever terms
    will be added if other datasets are incorporated.

    ---------------------------------------------------------------------------------------------------------------------

    Additional Fields:
    ------------------
    unique_topics: Unique discovered topics as determined within 'determine_topics' function.

    seed_terms: Optional BERTopic model seed terms (see bertopic package). Should be a list of dictionaries with {'Gene' : [term]} structure 
    or a list of terms. If a dictionary is supplied, the dictionary with a key matching the gene_name is selected.

    gene_text: Text linked to the gene which will be used to discover topics. If 'isolate_gene_text' function is not used, default to isolation 
    of 'Abstract' text in 'determine_topics' function. See 'isolate_gene_text' function for additional details.

    linked_topics: Topics which overlap between discovered topics and supplied topics of interest. See 'determine_topics' function.

    linked_genes: Gene names which are determined to be linked to the discovered topics. See 'determine_topics' function.
    '''
    def __init__(self, 
                 gene_name, 
                 description = None, 
                 pmid_data = ['PMID', 'Title', 'Abstract', 'Keywords', 'MeSH Terms']):
        
        self.gene_name = gene_name # Name of gene
        self.description = description # Desired gene description
        self.bertopic_model = None
        self.raw_topics = None
        self.unique_topics = None
        self.gene_text = []
        self.linked_genes = []
        self.topic_links = []
        self.pmid_data = {}

        for field in pmid_data:
            self.pmid_data[field] = [] 

    @classmethod           
    def preprocess_text(self, text):
        stop_words = stopwords.words('english')
        bound_terms = [r"\b" + term + r"\b" for term in stop_words]
        regex_pattern = "|".join(bound_terms)
        for i, string in enumerate(text):
            string = str(string)
            string = re.sub(r'\.(?=[a-zA-Z])', '. ', string) # If space not found between sentences, add one
            string = re.sub(r'\S*@\S*\s?', '', string)  # Remove emails
            string = re.sub('\'', '', string)  # Remove apostrophes
            string = re.sub(r'[^a-zA-Z0-9. ]+', ' ', string)  # Remove non-alphanumeric, periods, or spaces
            string = string.lower()  # Convert to lowercase
            string = re.sub(regex_pattern, '', string) # Remove stop words
            string = re.sub(r'\.[^.]*$', '', string) # Remove text after last period. This removes attribution for training purposes.
            string = re.sub(r'\s+', ' ', string)  # Remove extra spaces
            string = string.strip() # Remove trailing spaces
            text[i] = string
        return [item for item in text if item != ""]
        
    def __repr__(self):
        return f"{self.gene_name}"

    def __str__(self):
        return f"A Single_gene object named: {self.gene_name}"
        
    def add_PMID_data(self, pmid, search_terms = None, show_mismatch = False, check_en = False):
        '''
        pmid: List of dictionaries containing text data related to PMIDs of interest. Default data collected includes:
        PMID, Title, Keywords, and MeSH Terms.

        search_terms: List of search terms to use to ensure the linked PMID actually contains text related to the current gene.
        Defaults to None.

        show_mismatch: Bool, default to False. If set to True, compare text data extracted from current PMID to search terms. A 
        mismatch is defined as no search terms appearing in the PMID abstract text.

        check_en: Bool, default to False. If set to True, use langdetect to check each text section for text language.
        langdetect is a relatively computationally expensive process, so it is not enabled by default.
        '''

        if self.pmid_data.keys() != pmid.keys():
            missing_key = list(self.pmid_data.keys() - pmid.keys())
            for element in missing_key:
                pmid[element] = ''

        # Check if all of the pmid sections are in english
        if check_en == True:
            for key, value in pmid.items():
                if type(pmid[key]) == str and pmid[key] != '':
                    language = detect(pmid[key])
                    if language != 'en':
                        if show_mismatch == True:
                            print(f'Text was not english for {pmid['PMID']}. Skipping.')
                        return

        # Check if abstract text was found
        if pmid['Abstract'] == '':
             if show_mismatch == True:
                print(f'No abstract found for {pmid['PMID']}. Skipping.')
             return

        # If search terms are enabled, check if any of the search terms are in the found text
        if search_terms != None:
            if not isinstance(search_terms, list):
                    raise TypeError(f"search_terms expects a list of search terms. {type(search_terms)} was supplied.")
            else:
                pattern = r"\b(" + "|".join(re.escape(word) for word in search_terms) + r")\b"
                abstract_text = Single_gene.preprocess_text([pmid['Abstract']])
                if abstract_text == []: # For situations where no text is found after preprocessing, such as when abstract was not english.
                    return
                elif re.search(pattern, abstract_text[0]) == None:
                    if show_mismatch == True:
                        print(f'No matching terms found for {pmid['PMID']}. Skipping.')
                    return

        for key, value in pmid.items():
            if key not in self.pmid_data.keys():
                raise KeyError(f"{key} is not a valid data type. Currently available data types are: {list(self.pmid_data.keys())}")
            self.pmid_data[key].append(value)
    
    def set_bertopic(self, model, set_seed = False, seed_term_list = None):

        if set_seed == True:
            if isinstance(seed_term_list, list):
                seed_terms = seed_term_list
            elif seed_term_list == None:
                print('No seed_term_list was supplied. Defaulting to gene name.')
                seed_terms = self.gene_name.lower()
            else:
                raise TypeError(f'Supplied seed_term_list was {type(seed_term_list)}. Argument seed_term_list allows a list of search terms.')
            
            self.bertopic_model = model
            self.bertopic_model.seed_topic_list = seed_terms

    def isolate_gene_text(self, text_type = ['Abstract'], search_terms = None, debug = False):

        text_type_length = {k: len(v) for k, v in self.pmid_data.items() for k in text_type}
        if len(set(list(text_type_length.values()))) > 1:
            raise ValueError(f"Length of each self.pmid_data[text_type] should match but {len(set(list(text_type_length.values())))} lengths supplied. Check for dropped data from self.pmid.")
        else:
            text_length = max(set(list(text_type_length.values())))

        if self.description != None:
            self.gene_text.append(self.description[0]) # Add gene description to gene text

        if search_terms != None and not isinstance(search_terms, list):
            raise TypeError(f"search_terms expects a list of search terms. {type(search_terms)} was supplied.")

        self.gene_text = [] # Reset gene text stored in self
        found_gene_text = [] # Create storage list for gene text data

        # For each single paper stored in self.pmid_data, isolate and combine the text from that index
        for single_paper_data in range(text_length):

            all_text = [] # List to store extracted text

            for single_text_type in text_type:
                
                # From all pmid_data, isolate current text_type (i.e. Abstract text) from single_paper_data
                found_text = self.pmid_data[single_text_type][single_paper_data]

                if len(found_text) > 0:
                    found_text = found_text.strip() # Remove leading and trailing white space, if applicable.
                    if not found_text.endswith('.'): # If section does not end in a period, add a period
                        found_text += '.'
                    all_text.append(found_text)     
            
            if all_text != []:
                # Join text and check for a situation where space not added at start of sentence
                merged_text = ' '.join(all_text)
                merged_text = re.sub(r'\.(?=[a-zA-Z])', '. ', merged_text) 
                text = [item + '.' if not item.endswith('.') else item for item in merged_text.split('. ')] # Add period at end of sentences without a period.
                
                # If search_terms is not none, split the current text into individual 
                # sentences and retain only sentences containing search terms
                if search_terms != None:
                    sentences = [item for item in text if any(word in str(Single_gene.preprocess_text([item])) for word in search_terms)]
                    sentences = [x for x in sentences if re.search(r'[^\x00-\x7F]', x) == None] # Remove sentences with non-english characters
                    
                else:
                    sentences = [x for x in text if re.search(r'[^\x00-\x7F]', x) == None] # Remove sentences with non-english characters

                if sentences != []:
                    found_gene_text.append(' '.join(sentences))

                    if debug == True:
                        print(f'{self.pmid_data['PMID'][single_paper_data]}: {' '.join(sentences)}')

        self.gene_text = found_gene_text
                    
        return self.gene_text

    def determine_topics(self, gene_dict, strip_number_rear = 2):
        '''
        Determine topics appearing in text data related to a single gene. 

        Parameters:
        -----------
        gene_dict: Dictionary of {Gene:[Term]} structure with genes as keys and related terms, gene names, and/or gene products as 
        values.
        
        strip_number_rear: Number of characters to strip from the end of a discovered topic when comparing to supplied search terms. This
        allows for the identification of genes or gene products which may be designated as subunits, i.e. XXXa1, XXXp1. Default to 2.
        '''
        if self.bertopic_model == None:
            raise TypeError(f"No BERTopic model was supplied. Set a BERTopic model using set_bertopic().")
        
        if len(self.gene_text) == 0:
            print('No gene text was present. Defaulting to Abstract text. See isolate_gene_text() to change text type.')
            self.gene_text = self.isolate_gene_text(text_type = ['Abstract'])
            
        if not isinstance(gene_dict, dict):
            raise TypeError("gene_dict expects a dictionary with {'Gene' : [term]} structure.")

        search_terms = [values for subdict in gene_dict.values() for values in subdict]

        ############################

        local_gene_text = self.gene_text
        self.bertopic_model.fit_transform(local_gene_text)
        top_topics = self.bertopic_model.get_topic_info()

        embeds = self.bertopic_model.topic_embeddings_
        top_topics['Embeddings'] = [json.dumps(col.tolist()) for col in embeds]
        topic_distr, _ = self.bertopic_model.approximate_distribution(local_gene_text)
        top_topics['Topic_dist'] = [json.dumps(col.tolist()) for col in topic_distr.T]
        self.raw_topics = top_topics[['Topic', 'Count', 'Representation', 'Embeddings', 'Topic_dist']]
        
        top_topics = top_topics[top_topics.Topic != -1] # Remove outlier topics if using HDBSCAN
        topic_data = [item for sublist in top_topics['Representation'] for item in sublist]
        topic_data = self.preprocess_text(topic_data)
        unique_words = list(set(topic_data))

        ############################

        # Returns linked keys (gene names) and the overlap between found matches and original search terms
        linked_keys, overlapping_terms = find_genes(unique_words, 
                                                    gene_dict = gene_dict, 
                                                    search_terms = search_terms, 
                                                    strip_number_rear = strip_number_rear,
                                                    preprocessed = True)

        # Save found gene names and the overlap between matches and original search terms
        self.linked_genes = np.unique([key for sublist in linked_keys for key in sublist if key != self.gene_name]).tolist()
        self.topic_links = overlapping_terms

        # Isolate the found topics (item at position 1 of overlapping_terms) and remove them from the unique topics
        found_topics = Single_gene.preprocess_text([item[1] for item in overlapping_terms])
        self.unique_topics = list(set(unique_words) - set(found_topics))
        
        return self.unique_topics, self.linked_genes
