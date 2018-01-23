"""
Grand Bazar
"""
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from Bio import Entrez
from Bio.Entrez import efetch, read
import nltk
import re
import pprint
import glob
import os

from bioservices import *

import mygene


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


def save_abstract(abstract, save_file):
	##
	## -> Save the abstract in a text file
	## convert the abstract to unicode.
	##

	## preprocess abstract
	#abstract_preprocess = unicode(abstract)
	abstract_preprocess = abstract.encode('utf8')

	## save abstract in file
	output = open(save_file, "w")
	output.write(abstract_preprocess)
	output.close()

def get_cohorte_size(abstract_file):
	##
	## IN PROGRESS
	##
	## -> Try to retrieve the number of patients/case
	## in the dataset presented in the abstract
	## 
	## TODO : add more regex for size detection
	## TODO : add to bibotlite.py
	##

	## read abstract
	abstract = open(abstract_file, "r")

	cohorte_size = -1
	for line in abstract:

		## Try to catch a number before the apparition of the 
		## term in catch terms list. keep the biggest number
		## found in the abstract.
		catch_terms = ["patients", "cases", "observations", "subjects"]
		for item in catch_terms:
			m = re.findall(r"([0-9]+[\.|,]?[0-9]+) "+str(item), line)		
			if(m is not None):
				for match in m:
					match = match.replace(",", "")
					try:
						fetched_cohorte_size = int(match)
						if(fetched_cohorte_size > cohorte_size):
							cohorte_size = fetched_cohorte_size
					except:
						## do nothing
						tardis = 1

		## Try to catch a number after the apparition of the 
		## term "n = ". keep the biggest number found in the
		## abstract.
		m = re.findall(r"n = ([0-9]+[\.|,]?[0-9]+)", line)		
		if(m is not None):
			for match in m:
				match = match.replace(",", "")
				try:
					fetched_cohorte_size = int(match)
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
				except:
					## do nothing
					tardis = 1

		## Try to catch a number in a classic expression
		m = re.findall(r"[I,i]n total, ([0-9]+[\.|,]?[0-9]+) subjects", line)
		if(m is not None):
			for match in m:
				match = match.replace(",", "")
				try:
					fetched_cohorte_size = int(match)
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
				except:
					## do nothing
					tardis = 1

		"""
		## Try to catch a number in a classic expression
		m = re.findall(r"([0-9]+[\.|,]?[0-9]+)", line)
		if(m is not None):
			for match in m:
				match = match.replace(",", "")
				try:
					fetched_cohorte_size = int(match)
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
				except:
					## do nothing
					tardis = 1
		"""


	## close abstract
	abstract.close()

	## check if the cohorte_size variable
	## is still set to -1, set it to NA if
	## True
	if(cohorte_size == -1):
		cohorte_size = "NA"

	## return the detected size of the cohorte
	return cohorte_size





def get_date_from_meta_save(meta_file):
	##
	## Get the date of an article using the
	## meta data file created on local device,
	## no connection needed to NCBI server
	##
	## -> return the year of publication
	##
	## TODO : add to bibotlite.py
	##

	## Retrieve the year of publication
	## from the meta data file.
	year = "NA"
	meta_data = open(meta_file, "r")
	for line in meta_data:
		
		line = line.replace("\n", "")
		if(line[0] == ">"):

			line_in_array = line.split(";")
			if(line_in_array[0] == ">Date"):
				date_in_array = line_in_array[1].split("/")
				year = date_in_array[2]

	meta_data.close()

	## return only the year of publication
	return year




def plot_pulbications_years(meta_data_folder):
	##
	## Retrieve the year of publications of all
	## articles from the meta_data_folder and
	## plot the histogramm of publications over
	## the years
	##
	## TODO : add to bibotlite.py
	##

	## create the structure
	year_to_count = {}
	for meta_file in glob.glob(meta_data_folder+"/*.csv"):
		year = get_date_from_meta_save(meta_file)
		
		if(int(year) < 2018):

			if(year not in year_to_count.keys()):
				year_to_count[year] = 1
			else:
				year_to_count[year] += 1

	## plot graphe
	plt.bar(year_to_count.keys(), year_to_count.values(), 1, color='b')
	plt.show()




def plot_word_evolution(item_list, run_folder):
	##
	## -> Plot word frequency evolution over
	## the last decade
	##


	for item in item_list:

		## get year to pmid
		year_to_pmid = {}
		for meta_file in glob.glob(run_folder+"/meta/*.csv"):
			pmid = meta_file.split("/")
			pmid = pmid[-1]
			pmid = pmid.split(".")
			pmid = pmid[0]
			year = get_date_from_meta_save(meta_file)
			if(year not in year_to_pmid.keys()):
				year_to_pmid[year] = []
				year_to_pmid[year].append(pmid)
			else:
				year_to_pmid[year].append(pmid)

		## get apparition count in articles
		year_to_count = {}
		for year in year_to_pmid.keys():
			list_of_pmis_to_check = year_to_pmid[year]
			year_to_count[year] = 0
			for pmid in list_of_pmis_to_check:
				abstract_file = run_folder+"/"+"abstract/"+str(pmid)+"_abstract.txt"
				try:
					abstract_data = open(abstract_file, "r")
					for line in abstract_data:	
						m = re.findall(r"("+str(item)+")", line)		
						if(m is not None):
							if(len(m) > 0):
								year_to_count[year] += 1
					abstract_data.close()
				except:
					## do nothing
					tardis = 1

		## get apparition frequencies
		year_to_frequency = {}
		for year in year_to_pmid.keys():
			frequency = float(year_to_count[year])/float(len(year_to_pmid[year]))
			year_to_frequency[year] = frequency

		## save the figure
		y_vector = []
		x_vector = sorted(year_to_frequency.keys())

		for year in x_vector:
			y_vector.append(year_to_frequency[year])
		fig = plt.figure()
		plt.plot(x_vector, y_vector)
		plt.savefig(str(item)+"_frequency.png")
		#plt.show()
	




def describe_articles_type(run_folder):
	##
	## 
	## Count the number of articles talking about
	## several diseases, return a dictionnary.
	## currently work on:
	##		- SjS
	##		- SLE
	##		- RA
	##

	abstract_to_disease = {}
	abstract_to_disease["SjS"] = 0
	abstract_to_disease["SLE"] = 0
	abstract_to_disease["RA"] = 0 
	abstract_to_disease["Other"] = 0
	abstract_list = glob.glob(str(run_folder)+"/abstract/*.txt")
	for abstract_file in abstract_list:

		talking_about_SjS = False
		talking_about_RA = False
		talking_about_SLE = False

		abstract = open(abstract_file, "r")

		for line in abstract:

			## SjS
			m = re.findall(r"([S,s]j.gren)", line)
			if(len(m)>0):
				talking_about_SjS = True
			m = re.findall(r"(SjS)", line)
			if(len(m)>0):
				talking_about_SjS = True

			## SLE
			m = re.findall(r"([L,l]upus)", line)
			if(len(m)>0):
				talking_about_SLE = True
			m = re.findall(r"(SLE)", line)
			if(len(m)>0):
				talking_about_SLE = True

			## RA
			m = re.findall(r"([A,a])rthrisis", line)
			if(len(m)>0):
				talking_about_RA = True
			m = re.findall(r"(RA)", line)
			if(len(m)>0):
				talking_about_RA = True

		
		if(talking_about_SjS):
			abstract_to_disease["SjS"] += 1
		if(talking_about_SLE):
			abstract_to_disease["SLE"] += 1
		if(talking_about_RA):
			abstract_to_disease["RA"] += 1
		else:
			abstract_to_disease["Other"] += 1

			
		abstract.close()

	

	print abstract_to_disease
		







"""
TEST SPACE
"""
"""
print "------[TEST SPACE]------\n"
machin = get_ListOfDirectInteraction("P43403", "P06239")
print "-------------------------------------------------------"
print machin
"""

#plot_pulbications_years("SAVE/run_2/meta")


describe_articles_type("SAVE/run_2")

#item_list = ["neural network", "machine learning", "machine", "classification", "modelisation", "Sjogren", "random forest", "kmean",
#"statistic", "bioinformatic", "big data", "artificial intelligence", "diagnostic", "patients", "learning", "prediction", "cluster", "clusterring",
#"computer", "lupus", "RA", "IA", "sjogren"]
#plot_word_evolution(item_list, "SAVE/run_2")

"""
m = re.search(r"(?P<number>\w+) patients", "There are 257 patients in this study.")
if(m is not None):
	print m.group('number')
"""

"""
pmid_list = []
pmid_file = open("/home/perceval/Workspace/publications/immuno_review/articles/manually_retrieved/MANIFEST.txt", "r")
for line in pmid_file:
	line = line.replace("\n", "")
	pmid_list.append(line)
pmid_file.close()
for pmid in pmid_list:
	try:
		abstract = fetch_abstract(pmid)
		save_abstract(abstract, "abstract_1/abstract/"+str(pmid)+".txt")
	except:
		print "FAILED"
"""

"""
save_path = "SAVE/run_1/abstract/*.txt"
old_path = "abstract_1/abstract/*.txt"
run_2 = "abstract_2/*.txt"
current_run = "abstract/*.txt"
cmpt = 1
for abstract in glob.glob(current_run):
	size = get_cohorte_size(abstract)
	print size
	cmpt += 1
"""


#k = KEGG(verbose=False)
#truc = k.find("hsa", "zap70")
#pathway = k.get_pathway_by_gene("1956", "hsa")
#print pathway
#k.show_pathway("hsa05218", keggid={"1956": "red"})



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
"""
print "[+] => Testing get_ListOfArticles function"
machin = get_ListOfArticles("HLA", 100)
for pmid in machin:
	try:
		test = fetch_abstract(pmid)
		print "["+str(pmid)+"]\n"
		print test
	except:
		print "[!] Can't read "+str(pmid)
print "[*] => Test Done"
"""


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

## Test mygene module
"""
mg = mygene.MyGeneInfo()
truc = mg.query('NRAS', size=1)
print truc["hits"][0]
"""