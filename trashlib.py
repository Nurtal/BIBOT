"""
Grand Bazar
"""

from Bio import Entrez
from Bio.Entrez import efetch, read
import nltk
import re
import pprint

import os

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
	xml_data = read(handle)
	xml_data = xml_data['PubmedArticle'][0]
	
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



def get_interactors(interaction):
	"""
	return the list of elements in interaction
	-> always 2 elements
	-> work for mint db
	"""

	list_ofInteractors = []

	line1 = interaction[0]
	line2 = interaction[1]
	line1InArray = line1.split(":")
	line2InArray = line2.split(":")
	if(len(line1InArray) > 1):
		uniprotId_elt1 = line1InArray[1]
		list_ofInteractors.append(str(uniprotId_elt1))
	else:
		uniprotId_elt1 = "CAN'T PARSE DATA"
		list_ofInteractors.append(str(uniprotId_elt1))
	if(len(line2InArray) > 1):
		uniprotId_elt2 = line2InArray[1]
		list_ofInteractors.append(str(uniprotId_elt2))
	else:
		uniprotId_elt2 = "CAN'T PARSE DATA"
		list_ofInteractors.append(str(uniprotId_elt2))

	return list_ofInteractors



def get_CuratedInteraction(elementToCheck):
	"""
	return the list of curated interaction (i.e doublons are removed)
	including elementToCheck.

	-> work for mint db
	-> id are uniprotId
	"""
	#list_elementsToCheck.append(str(elementToCheck))
	s = PSICQUIC(verbose=False)
	data = s.query("mint", str(elementToCheck)+" AND species:9606")

	listOfInterctorsInGraph = []
	listOfCuratedInteraction = []

	for db in data:
		machin = get_interactors(db)
		interactionType = get_InteractionType(db)
		if(str(machin[0]) == str(elementToCheck) and str(machin[1]) not in listOfInterctorsInGraph):
			interaction = str(machin[0])+"\t"+str(interactionType)+"\t"+str(machin[1])
			listOfCuratedInteraction.append(interaction)
			listOfInterctorsInGraph.append(machin[1])
		elif(str(machin[1]) == str(elementToCheck) and str(machin[0]) not in listOfInterctorsInGraph):
			interaction = str(machin[1])+"\t"+str(interactionType)+"\t"+str(machin[0])
			listOfCuratedInteraction.append(interaction)
			listOfInterctorsInGraph.append(machin[0])

	return listOfCuratedInteraction




def draw_InteractionGraph(element, fileName):
	"""
	Write all interactions including element and interactions
	including each element assiocated with the first element
	and so on ... in a .sif file (fileName).

	-> Work with mint database
	-> use grep -v "sapiens" to deal with unwanted output
	return nothing
	"""

	list_elementsToCheck = []
	list_elementChecked = []
	list_elementsToCheck.append(element)

	for element in list_elementsToCheck:
		if(element not in list_elementChecked):
			machin = get_CuratedInteraction(str(element))
		
			list_elementsToCheck.remove(str(element))
			list_elementChecked.append(str(element))

			for interaction in machin:
				interactionInArray = interaction.split("\t")

				grapheFile = open(fileName, "a")
				grapheFile.write(interaction+"\n")
				grapheFile.close()

				print interaction
				if(interactionInArray[2] not in list_elementsToCheck):
					list_elementsToCheck.append(interactionInArray[2])
		else:
			print "->"+str(element) + "already Ckeck"



def get_InteractionType(interaction):
	"""
	return the type of interaction
	-> work on mint db
	"""

	typeOfInteraction =  interaction[11]
	typeOfInteractionInArray = typeOfInteraction.split("(")
	interactionName = typeOfInteractionInArray[1]
	interactionName = interactionName[:-1]
	return interactionName



def convert_SifFileToGDFfile(fileName):
	"""
	[IN PROGRESS]

	-> Try on small files 

	"""

	#check extension
	fileNameInArray = fileName.split(".")
	if(len(fileNameInArray) < 1):
		print "Can't parse format for conversion\n"
	elif(fileNameInArray[1] != "sif"):
		print "Wrong format for input file, .sif expected ("+fileNameInArray[1]+" found)"
	else:
		inputFilename = fileNameInArray[0]

	#Initialise node tmp file
	tmpFile_node = open("nodes_buildGDF.tmp", "a")
	tmpFile_node.write("nodedef>name VARCHAR\n")
	tmpFile_node.close()

	#Initialise edge tmp file
	tmpFile_edge = open("edges_buildGDF.tmp", "a")
	tmpFile_edge.write("edgedef>node1 VARCHAR,node2 VARCHAR\n")
	tmpFile_edge.close()

	listOfNodes = []
	sifFile = open(fileName, "r")
	for line in sifFile:
		lineInArray = line.split("\t")
		node1 = lineInArray[0]
		node2 = lineInArray[2]
		if(node1 not in listOfNodes):
			tmpFile_node = open("nodes_buildGDF.tmp", "a")
			tmpFile_node.write(str(node1)+"\n")
			tmpFile_node.close()
			listOfNodes.append(node1)
		
		tmpFile_edge = open("edges_buildGDF.tmp", "a")
		tmpFile_edge.write(str(node1)+","+str(node2))
		tmpFile_edge.close()

	# Write GDF file
	outputFile = open(str(inputFilename)+".gdf", "a")
	tmpFile_node = open("nodes_buildGDF.tmp", "r")
	for line in tmpFile_node:
		outputFile.write(line)
	tmpFile_node.close()
	tmpFile_edge = open("edges_buildGDF.tmp", "r")
	for line in tmpFile_edge:
		outputFile.write(line)
	tmpFile_edge.close()
	outputFile.close()

	# delete tmp file
	os.remove("nodes_buildGDF.tmp")
	os.remove("edges_buildGDF.tmp")





"""
TEST SPACE
"""
print "------[TEST SPACE]------\n"

"""
machin = get_ListOfDirectInteraction("P43403", "P06239")
print "-------------------------------------------------------"
print machin
"""


"""
k = KEGG(verbose=False)
k.find("hsa", "zap70")
pathway = k.get_pathway_by_gene("7535", "hsa")
print pathway
k.show_pathway("hsa04064", keggid={"7535": "red"})
"""



#draw_InteractionGraph("P43403", "monTest.sif")
#convert_SifFileToGDFfile("monTest.sif")

#abstract = fetch_abstract(27755966)
#print abstract



"""
list_elementsToCheck = []
list_elementChecked = []
list_elementsToCheck.append("P43403")

for element in list_elementsToCheck:
	if(element not in list_elementChecked):
		machin = get_CuratedInteraction(str(element))
		
		list_elementsToCheck.remove(str(element))
		list_elementChecked.append(str(element))

		for interaction in machin:
			interactionInArray = interaction.split("->")
			print interaction
			if(interactionInArray[1] not in list_elementsToCheck):
				list_elementsToCheck.append(interactionInArray[1])
	else:
		print "->"+str(element) + "already Ckeck"
"""



"""
elementToCheck = "ZAP70"
s = PSICQUIC(verbose=False)
#data = s.query("mint", "ZAP70 AND species:9606")
data = s.query("mint", str(elementToCheck)+" AND species:9606")

for db in data:
	machin = get_interactors(db)
	truc = get_InteractionType(db)
	print truc
"""

	

	

"""
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


## Test get_ListOfArticles function
print "[+] => Testing get_ListOfArticles function"
machin = get_ListOfArticles("Lymphoma", 1000)
for pmid in machin:
	try:
		test = fetch_abstract(pmid)
		print "["+str(pmid)+"]\n"
		print test
	except:
		print "[!] Can't read "+str(pmid)
print "[*] => Test Done"



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
