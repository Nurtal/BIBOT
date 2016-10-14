"""
Grand Bazar
"""

from Bio import Entrez
from Bio.Entrez import efetch, read
import nltk
import re
import pprint

from bioservices import *

Entrez.email = 'murlock.raspberypi@gmail.com'


def fetch_abstract(pmid):
	"""
	Retrun abstract of a given
	article using pmid

	=> Return None when pmid can't be return
	(can happen when article is in chinese)
	"""
	handle = efetch(db='pubmed', id=pmid, retmode='xml', )
	xml_data = read(handle)[0]

	try:
		article = xml_data['MedlineCitation']['Article']
		abstract = article['Abstract']['AbstractText'][0]
		return abstract
	except IndexError:
		return None
	except KeyError:
		return None


def fetch_article(pmid):
	"""
	Test function
	"""
	handle = efetch(db='pubmed', id=pmid, retmode='xml', )
	xml_data = read(handle)[0]

	try:
		article = xml_data['MedlineCitation']['Article']
		abstract = article['Abstract']['AbstractText'][0]
		return article

	except IndexError:
		return None


def get_ListOfArticles(term, max_count):
	"""
	return the list of pmid article conatining
	the term.
	"""
	h = Entrez.esearch(db='pubmed', retmax=max_count, term=term)
	result = Entrez.read(h)
	listOfArticles = result["IdList"]


	return listOfArticles;



def get_ListOfCommunPathway(elt1, elt2):

	"""
	-> elt1 and elt2 are KEGG id
	return the list of commun pathway between elt1 and elt2
	-> Use bioservices
	-> Internet connection needed
	-> Use KEGG database
	-> "hsa" correspond to human in KEGG database
	"""

	k = KEGG(verbose=False)

	list_OfCommunPathway = []
	list_pathways_elt1 = k.get_pathway_by_gene(str(elt1), "hsa")
	list_pathways_elt2 = k.get_pathway_by_gene(str(elt2), "hsa")

	for pathway_elt1 in list_pathways_elt1.keys():
		for pathway_elt2 in list_pathways_elt2.keys():
			if(str(pathway_elt1) == str(pathway_elt2)):
				list_OfCommunPathway.append(str(pathway_elt1))

	return list_OfCommunPathway



def get_ListOfDirectInteraction(elta, eltb):
	"""
	-> Scan interaction DataBase and return the list of direct interaction
	between elta and eltb.
	-> For now elta and eltb have to be uniprot ID
	-> For now only scan mint DB
	"""

	list_interactions = []
	s = PSICQUIC(verbose=False)
	data = s.query("mint", str(elta)+" AND species:9606")

	for interaction in data:
		line1 = interaction[0]
		line2 = interaction[1]
		line1InArray = line1.split(":")
		line2InArray = line2.split(":")
		uniprotId_elt1 = line1InArray[1]
		uniprotId_elt2 = line2InArray[1]
		if(str(uniprotId_elt1) == str(elta) and str(uniprotId_elt2) == str(eltb)):
			list_interactions.append(interaction)


	return list_interactions







"""
TEST SPACE
"""
print "------[TEST SPACE]------\n"


machin = get_ListOfDirectInteraction("P43403", "P06239")
print "-------------------------------------------------------"
print machin



"""
k = KEGG(verbose=False)
k.find("hsa", "zap70")
pathway = k.get_pathway_by_gene("7535", "hsa")
print pathway
k.show_pathway("hsa04064", keggid={"7535": "red"})
"""

"""
s = PSICQUIC(verbose=False)
#data = s.query("mint", "ZAP70 AND species:9606")
data = s.query("mint", "P43403 AND species:9606")

for db in data:
	line1 = db[0]
	line2 = db[1]
	line1InArray = line1.split(":")
	line2InArray = line2.split(":")
	uniprotId_elt1 = line1InArray[1]
	uniprotId_elt2 = line2InArray[1]
	print str(uniprotId_elt1)+" || "+str(uniprotId_elt2)
"""





# Retrieve Article
#test = fetch_abstract(27045581)
#print test

#machin = get_ListOfArticles("B-Cell", 1000)
#for pmid in machin:
#	test = fetch_abstract(pmid)
#	print "["+str(pmid)+"]\n"
#	print test




"""
# test data
test = "Rabbits are dangerous. Rabbits are not dangerous"
print test
# import nltk stuff
nltk.download('punkt')
nltk.download('averaged_perceptron_tagger')

# Structure Information
sentences = nltk.sent_tokenize(test)
sentences = [nltk.word_tokenize(sent) for sent in sentences]
sentences = [nltk.pos_tag(sent) for sent in sentences]

# Chunking
grammar = "NP: {<DT>?<JJ>*<NN>}"
cp = nltk.RegexpParser(grammar)
for sentence in sentences:
	print sentence
	result = cp.parse(sentence)
	result.draw()
"""
