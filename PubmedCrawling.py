# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:01:10 2024

@author: HONGYAN
"""

import csv
from Bio import Entrez

# Always tell NCBI who you are (your email address)
Entrez.email = "your.email@example.com"

def search_pubmed(query, max_results=30):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax=max_results,
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    handle.close()
    return results['IdList']

def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    handle.close()
    return results

def save_to_csv(articles, filename='pubmed_results.csv'):
    with open(filename, 'a', newline='', encoding='utf-8') as csvfile:  # 'a' for append mode
        fieldnames = ['PMID', 'Title', 'Authors', 'Journal', 'Abstract', 'PubDate']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Check if file is empty to write headers or not
        csvfile.seek(0, 2)  # Move to the end of the file
        if csvfile.tell() == 0:
            writer.writeheader()

        for article in articles['PubmedArticle']:
            article_data = article['MedlineCitation']['Article']
            pmid = article['MedlineCitation']['PMID']
            title = article_data.get('ArticleTitle', '')
            abstract = article_data.get('Abstract', {}).get('AbstractText', [''])[0]
            authors_list = article_data.get('AuthorList', [])
            authors = ', '.join([' '.join([author.get('LastName', ''), author.get('Initials', '')]) for author in authors_list])
            journal = article_data.get('Journal', {}).get('Title', '')
            pub_date = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', '')
            writer.writerow({'PMID': pmid, 'Title': title, 'Authors': authors, 'Journal': journal, 'Abstract': abstract, 'PubDate': pub_date})

if __name__ == '__main__':
    queries = ['cancer', 'sex difference']
    for query in queries:
        id_list = search_pubmed(query)
        articles = fetch_details(id_list)
        save_to_csv(articles)
    print(f"Results saved to 'pubmed_results.csv'")
