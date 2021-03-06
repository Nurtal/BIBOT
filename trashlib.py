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
import operator

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
		





def get_all_subjects_from_abstract(abstract_folder):
	##
	## The idea is to find all key word from
	## each abstract in a run folder and
	## return a list of key word.
	##
	## Word are selected by the number of apparition in
	## all abstracts, if found in more than 10 % of the articles,
	## the word is selected
	##

	## get list of all abstract files
	abstract_list = glob.glob(abstract_folder+"/*.txt")
	name_list = []
	name_to_count = {}

	treshold = 0.1

	abstract_parsed = 0
	for abstract in abstract_list:
		names_found_in_abstract = []
		abstract_data = open(abstract, "r")
		for line in abstract_data:
			try:
				tokens = nltk.word_tokenize(line.encode('utf8'))
				tagged = nltk.pos_tag(tokens)
				entities = nltk.chunk.ne_chunk(tagged)
				abstract_parsed += 1
			except:
				## Something went wrong
				entities = []

			for item in entities:
				try:
					if(item[1] in ["NN", "NNS", "NNP"]):
						if(item[0] not in names_found_in_abstract):
							names_found_in_abstract.append(item[0])
				except:
					## Somethig went wrong
					choucroute = True

		for name in names_found_in_abstract:
			if(name not in name_to_count.keys()):
				name_to_count[name] = 1
			else:
				name_to_count[name] += 1
		abstract_data.close()

	for name in name_to_count.keys():
		if(float(name_to_count[name]/float(abstract_parsed)) > treshold and name not in name_list):
			name_list.append(name)	

	return name_list



def get_number_of_articles_from_log_file(log_file):
	##
	## Get the total number of articles found
	## by the run by looking to the run's log
	## file.
	##
	## return an int or "NA" if nothing is found. 
	##

	number_of_articles = "NA"
	log_data = open(log_file, "r")
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(line_in_array[0] == "Total_number_of_articles"):
			try:
				number_of_articles = int(line_in_array[1])
			except:
				number_of_articles = "NA"

	log_data.close()

	return number_of_articles


def get_number_of_keywords_from_log_file(log_file):
	##
	## Get the total number of keywords used to
	## generate the combination of request by
	## the run.
	##
	## return an int or "NA" if nothing is found. 
	##

	number_of_keywords = "NA"
	log_data = open(log_file, "r")
	cmpt = 0
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(cmpt == 0):
			try:
				number_of_keywords = len(line_in_array)
			except:
				number_of_keywords = "NA"
		else:
			break

		cmpt += 1

	log_data.close()

	return number_of_keywords



def get_articles_filter_stat(log_file):
	##
	## Get the count of articles that passed the first
	## filter, the second filter, and both.
	## return a tuple of int
	##

	first_filter_pass_cmpt = 0
	second_filter_pass_cmpt = 0
	pass_both_filter_cmpt = 0

	log_data = open(log_file, "r")
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(line[0] == ">"):
			
			filter_1_info = line_in_array[1].split("=")
			if(filter_1_info[1] == "PASSED"):
				first_filter_pass_cmpt += 1

			filter_2_info = line_in_array[2].split("=")
			if(filter_2_info[1] == "PASSED"):
				second_filter_pass_cmpt += 1


			if(filter_1_info[1] == "PASSED" and filter_2_info[1] == "PASSED"):
				pass_both_filter_cmpt += 1


	log_data.close()

	return (first_filter_pass_cmpt, second_filter_pass_cmpt, pass_both_filter_cmpt)



def get_number_of_articles_for_year(run_folder, year):
	##
	## get the number of articles published in a
	## specific year, parse data from meta files
	## in the given run folder
	##

	article_cmpt = 0

	for meta_file in glob.glob(run_folder+"/meta/*.csv"):
		meta_data = open(meta_file, "r")
		for line in meta_data:
			line = line.replace("\n", "")
			line_in_array = line.split(";")
			if(line_in_array[0] == ">Date"):
				date = line_in_array[1]
				date_in_array = date.split("/")
				year_fetched = date_in_array[-1]
				if(int(year_fetched) == year):
					article_cmpt += 1
		meta_data.close()

	return article_cmpt


def get_country_publication_stat(run_folder):
	##
	## Get the publications stats for country.
	## get the list of pmid retrieved from the
	## meta folder and connect to the NCBI to fecth
	## publications informations, parse it to get the
	## country of publication.
	## 
	## return a dictionnary
	##

	## init structure
	country_to_count = {}

	## get list of PMID to process
	meta_file_list = glob.glob(run_folder+"/meta/*.csv")
	for meta_file in meta_file_list:
		meta_file_in_array = meta_file.split("/")
		file_name = meta_file_in_array[-1]
		file_name_in_array = file_name.split(".")
		pmid = file_name_in_array[0]

		## get country publication
		try:
			handle = efetch(db='pubmed', id=pmid, retmode='xml', )
			informations = read(handle)
			stuff = informations[u'PubmedArticle'][0]
			country = stuff[u'MedlineCitation'][u'MedlineJournalInfo'][u'Country']
		except:
			country = "NA"

		## fill dictionnary
		if(country not in country_to_count.keys()):
			country_to_count[country] = 1
		else:
			country_to_count[country] += 1

	return country_to_count





def get_article_domain(abstract):
	##
	## IN PROGRESS
	##
	## Fit an article on a pre existing
	## class based on the content of its abstract
	##
	##
	## Class:
	##
	##	- Diagnostic : predict the diagnostic of a patient
	##
	##	- Therapeutic : test a treatment
	##
	##	- Modelisation : try to understand the disease
	##
	## 	- Unclassified : not sure what we are talking about
	##

	## Control structure
	diagnostic_keywords = ["classification", "Classification", "criteria", "diagnostic", "diagnosis", "prevalence", "epidemiological"]
	therapeutic_keywords = ["therapeutic", "therapy", "treatment", "treatments", "rituximab"]
	modelisation_keywords = ["model", "models", "modelisation", "components", "dynamics", "composition", "pathway", "regulatory", "regulates", "mechanistic", "mechanism", "mechanisms"]

	diagnostic_score = 0
	therapeutic_score = 0
	modelisation_score = 0

	article_theme = "Unclassified"

	## Looking for keyword in the abstract with nltk
	names_found_in_abstract = []
	abstract_data = open(abstract, "r")
	for line in abstract_data:
		try:
			tokens = nltk.word_tokenize(line.encode('utf8'))
			tagged = nltk.pos_tag(tokens)
			entities = nltk.chunk.ne_chunk(tagged)
		except:
			## Something went wrong
			entities = []
		for item in entities:
			try:
				if(item[1] in ["NN", "NNS", "NNP"]):
					if(item[0] not in names_found_in_abstract):
						names_found_in_abstract.append(item[0])
			except:
				## Somethig went wrong
				choucroute = True

	## compute score
	for name in names_found_in_abstract:
		if(name in diagnostic_keywords):
			diagnostic_score += 1
		elif(name in therapeutic_keywords):
			therapeutic_score += 1
		elif(name in modelisation_keywords):
			modelisation_score += 1

	abstract_data.close()
	
	## looking for regular expression
	## split abstract into sentences, look for combination
	## of keywords in each sentences.
	## Get the list of words
	abstract_data = open(abstract, "r")
	text = ""
	for line in abstract_data:
		text +=line
	abstract_data.close()
	sentences = text.split(". ")
	words = text.replace(".", "")
	words = words.replace(",", "")
	words = words.replace(":", "")
	words = words.replace(";", "")
	words = words.split(" ")


	## Looking for co-occurences
	for sentence in sentences:
		
		## modelisation 
		if(("mechanisms" in sentence or "mechanism" in sentence) and ("understood" in sentence or "understand" in sentence)):
			modelisation_score += 1
		if(("summarize" in sentence) and ("mechanism" in sentence or "mechanisms" in sentence or "pathway" in sentence or "pathways" in sentence)):
			modelisation_score += 1
		if(("present" in sentence and "study" in sentence) and "investigated" in sentence):
			modelisation_score += 1
		if("investigated" in sentence and "interactions" in sentence):
			modelisation_score += 1

		## Diagnostic
		if("identification" in sentence and ("phenotype" in sentence or "phenotypes" in sentence or "subphenotype" in sentence or "subphenotypes" in sentence)):
			diagnostic_score += 1
		if("severity" in sentence and "marker" in sentence):
			diagnostic_score += 1
		if("differentiate" in sentence and "patients" in sentence):
			diagnostic_score += 1
		if(("biomarker" in sentence or "biomarkers" in sentence) and ("disease" in sentence or "diseases" in sentence)):
			diagnostic_score += 1

		## Therapeutic
		if("impact" in sentence and "course of disease" in sentence):
			therapeutic_score += 1

	## Checking words
	for word in words:
		if(word in diagnostic_keywords):
			diagnostic_score += 1
		elif(word in therapeutic_keywords):
			therapeutic_score += 1
		elif(word in modelisation_keywords):
			modelisation_score += 1



	## compute score
	scores = {"Diagnostic":diagnostic_score, "Therapeutic":therapeutic_score, "Modelisation":modelisation_score}
	sorted_scores = sorted(scores.items(), key=operator.itemgetter(1))

	if(sorted_scores[-1][1] > 0):
		article_theme = sorted_scores[-1][0]

	return article_theme



def get_count_of_articles_type_for_each_disease():
	##
	## IN PROGRESS
	##

	## Diagnostic
	for file in glob.glob("articles_classification/Diagnostic/*"):

		abstract = open(file, "r")
		abstract_in_line = ""
		for line in abstract:
			abstract_in_line += line
		abstract.close()

		## SjS
		disease_keywords = ["SjS", "sjogren"]




def retrieve_techniques(display_details):
	##
	## Retrieve techniques used in articles from parsing the abstracts
	## display_details is a boolean:
	##    - True: display details
	##    - False: do nothing 
	##
	##
	## Read abstracts in the articles_classification subfolders
	##    - Diagnostic
	##    - Therapeutic
	##    - Modelisation
	##    - Unclassified
	## 
	## Parse each abstract and look for a list of statistical analaysis appraoches:
	##    - t-test
	##    - logistic regression
	##    - cox regression
	##    - multivariate regression
	##    - univariate regression
	##    - linear regression
	##    - Poisson regression
	##    - regression tree
	##    - regression analysis
	##    - clustering
	##    - hierarchical clustering
	##    - decision tree
	##    - random forest
	##    - artificial neural network
	##    - support vector machine
	##    - machine learning
	##    - Wilcoxon test
	##    - PCA
	##    - chi-squared
	##    - multivariate analysis
	##    - univariate analysis
	##    - discriminant analysis
	##    - bioinformatics analysis
	##    - GLM (generalized linear model)
	##    - mixed models
	##    - normalization
	##    - k mean
	##
	## Return a dictionnary category : technique : count
	##


	## initiate data structure
	category_to_techniques = {"Diagnostic":{}, "Therapeutic":{}, "Modelisation":{}, "Unclassified":{}}

	## Parse the abstracts
	for category in ["Diagnostic", "Therapeutic", "Modelisation", "Unclassified"]:
		for file in glob.glob("articles_classification/"+str(category)+"/*"):

			abstract = open(file, "r")
			abstract_in_line = ""
			for line in abstract:
				abstract_in_line += line
			abstract.close()

			abstract_in_array = abstract_in_line.split(" ")

			cmpt = 0
			for word in abstract_in_array:

				## t-test
				if(word == "t-test"):
					if("t-test" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["t-test"] = 1
					else:
						category_to_techniques[str(category)]["t-test"] += 1
				elif(word == "t" and cmpt + 1 < len(abstract_in_array)):
					if(abstract_in_array[cmpt+1] in ["test", "tests"]):
						if("t-test" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["t-test"] = 1
						else:
							category_to_techniques[str(category)]["t-test"] += 1

				## regression
				##    - logistic
				##    - cox
				##    - multivariate
				##    - univariate
				##    - linear
				##    - tree
				##    - Poisson
				if(word in ["regression", "Regression"] and cmpt-1 >= 0):

					## logistic
					if(abstract_in_array[cmpt-1] in ["logistic", "Logistic", "logistical", "Logistical"]):
						if("logistic regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["logistic regression"] = 1
						else:
							category_to_techniques[str(category)]["logistic regression"] += 1

					## multivariate
					elif(abstract_in_array[cmpt-1] in ["multivariate", "Multivariate", "multiple", "Multiple"]):
						if("multivariate regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["multivariate regression"] = 1
						else:
							category_to_techniques[str(category)]["multivariate regression"] += 1
					
					## univariate
					elif(abstract_in_array[cmpt-1] in ["simple", "Simple", "univariate", "Univariate"]):
						if("univariate regression" not in category_to_techniques[str(category)]):
							category_to_techniques[str(category)]["univariate regression"] = 1
						else:
							category_to_techniques[str(category)]["univariate regression"] += 1

					## linear
					elif(abstract_in_array[cmpt-1] in ["Linear", "linear"]):
						if("linear regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["linear regression"] = 1
						else:
							category_to_techniques[str(category)]["linear regression"] += 1

					## cox
					elif(abstract_in_array[cmpt-1] in ["Cox", "cox", "Cox's", "cox's"]):
						if("Cox regression" not in category_to_techniques["Diagnostic"].keys()):
							category_to_techniques[str(category)]["Cox regression"] = 1
						else:
							category_to_techniques[str(category)]["Cox regression"] += 1

					## proportion hazard regression model
					elif(abstract_in_array[cmpt-1] in ["hazard", "hazards"]):
						if(abstract_in_array[cmpt-2] in ["proportion", "proportional"] and abstract_in_array[cmpt+1] in ["analysis", "model"]):
							if("proportion hazard regression model" not in category_to_techniques[str(category)].keys()):
								category_to_techniques[str(category)]["proportion hazard regression model"] = 1
							else:
								category_to_techniques[str(category)]["proportion hazard regression model"] += 1
					
					## regression tree
					elif(abstract_in_array[cmpt+1] == "tree"):
						if("regression tree" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["regression tree"] = 1
						else:
							category_to_techniques[str(category)]["regression tree"] += 1
					
					## Poisson regression
					elif(abstract_in_array[cmpt-1] in ["Poisson", "poisson"]):
						if("Poisson regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["Poisson regression"] = 1
						else:
							category_to_techniques[str(category)]["Poisson regression"] += 1

					elif(abstract_in_array[cmpt+1] in ["analysis", "test", "tests"]):
						if("regression analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["regression analysis"] = 1
						else:
							category_to_techniques[str(category)]["regression analysis"] += 1

				## Cox's hazard model
				elif(word in ["Cox's", "cox's", "cox", "Cox"] and abstract_in_array[cmpt+1] in ["proportional", "hazards", "hazard", "model", "models"]):
					if("proportion hazard regression model" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["proportion hazard regression model"] = 1
					else:
						category_to_techniques[str(category)]["proportion hazard regression model"] += 1

				## clustering
				elif(word in ["clustering", "Clustering"]):
					if(abstract_in_array[cmpt-1] in ["Hierarchical", "hierarchical"]):
						if("hierarchical clustering" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["hierarchical clustering"] = 1
						else:
							category_to_techniques[str(category)]["hierarchical clustering"] += 1
					else:
						if("clustering" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["clustering"] = 1
						else:
							category_to_techniques[str(category)]["clustering"] += 1

				## tree
				elif(word in ["tree", "Tree"]):
					if(abstract_in_array[cmpt-1] in ["decision", "Decision"]):
						if("decision tree" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["decision tree"] = 1
						else:
							category_to_techniques[str(category)]["decision tree"] += 1

				## random forest
				elif(word in ["forest", "Forest"]):
					if(abstract_in_array[cmpt-1] in ["random", "Random"]):
						if("random forest" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["random forest"] = 1
						else:
							category_to_techniques[str(category)]["random forest"] += 1

				## artificial neural network
				elif(word in ["neural", "Neural"]):
					if(abstract_in_array[cmpt+1] in ["network", "networks"] and (abstract_in_array[cmpt-1] in ["artificial", "Artificial"] or abstract_in_array[cmpt+2] in ["algorithm", "algorithms"])):
						if("artificial neural network" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["artificial neural network"] = 1
						else:
							category_to_techniques[str(category)]["artificial neural network"] += 1

				## machine learning
				## support vector machine
				elif(word in ["machine", "Machine"]):
					if(abstract_in_array[cmpt+1] in ["learning", "Learning"]):
						if("machine learning" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["machine learning"] = 1
						else:
							category_to_techniques[str(category)]["machine learning"] += 1

					elif(abstract_in_array[cmpt-1] == "vector" and abstract_in_array[cmpt-2] in ["support", "Support"]):
						if("support vector machine" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["support vector machine"] = 1
						else:
							category_to_techniques[str(category)]["support vector machine"] += 1
				elif(word in ["machine-learning", "Machine-learning"]):
					if("machine learning" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["machine learning"] = 1
					else:
						category_to_techniques[str(category)]["machine learning"] += 1

				## Wilcoxon test
				elif(word in ["Wilcoxon", "wilcoxon"]):
					if("Wilcoxon test" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["Wilcoxon test"] = 1
					else:
						category_to_techniques[str(category)]["Wilcoxon test"] += 1

				## PCA
				elif(word in ["Principal", "principal"]):
					if(abstract_in_array[cmpt+1] == "component" and abstract_in_array[cmpt+2] == "analysis"):
						if("PCA" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["PCA"] = 1
						else:
							category_to_techniques[str(category)]["PCA"] += 1

				## chi-squared
				elif(word in ["chi-squared", "Chi-squared", "Chi-sqare", "ch-square"]):
					if("chi-squared" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["chi-squared"] = 1
					else:
						category_to_techniques[str(category)]["chi-squared"] += 1

				## Fisher's test
				elif(word in ["Fisher's", "Fisher"]):
					if("Fisher test" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["Fisher test"] = 1
					else:
						category_to_techniques[str(category)]["Fisher test"] += 1

				## analysis
				## - multivariate
				## - univariate
				## - discriminant
				## - bioinformatics
				elif(word in ["analysis", "Analysis"]):

					if(abstract_in_array[cmpt-1] in ["multivariate", "Multivariate"]):
						if("multivariate analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["multivariate analysis"] = 1
						else:
							category_to_techniques[str(category)]["multivariate analysis"] += 1
					elif(abstract_in_array[cmpt-1] in ["univariate", "Univariate"]):
						if("univariate analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["univariate analysis"] = 1
						else:
							category_to_techniques[str(category)]["univariate analysis"] += 1
					elif(abstract_in_array[cmpt-1] in ["discriminant", "Discriminant"]):
						if("discriminant analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["discriminant analysis"] = 1
						else:
							category_to_techniques[str(category)]["discriminant analysis"] += 1

					elif(abstract_in_array[cmpt-1] in ["bioinformatics", "Bioinformatics", "bioinformatic", "Bioinformatic"]):
						if("bioinformatics analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["bioinformatics analysis"] = 1
						else:
							category_to_techniques[str(category)]["bioinformatics analysis"] += 1

				## GLM (generalized linear model)
				elif(word in ["GLM", "(GLM)", "(GLM, ", "GLM)"]):
					if("GLM" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["GLM"] = 1
					else:
						category_to_techniques[str(category)]["GLM"] += 1

				## mixed models
				elif(word in ["models", "model"] and abstract_in_array[cmpt-1] in ["mixed", "Mixed"]):
					if("mixed models" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["mixed models"] = 1
					else:
						category_to_techniques[str(category)]["mixed models"] += 1

				## normalization
				elif(word in ["Normalization", "normalization"]):
					if("normalization" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["normalization"] = 1
					else:
						category_to_techniques[str(category)]["normalization"] += 1

				## k mean
				elif(word in ["kmean", "Kmean", "K-mean", "K-mean-clustering"]):
					if("k mean" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["k mean"] = 1
					else:
						category_to_techniques[str(category)]["k mean"] += 1
				elif(word in ["k", "K"] and abstract_in_array[cmpt+1] in ["mean", "means", "Mean", "Means"]):
					if("k mean" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["k mean"] = 1
					else:
						category_to_techniques[str(category)]["k mean"] += 1

				cmpt += 1

	## Display details of the returned structure
	if(display_details):
		total = 119 + 604 +726 +19
		detected = 0
		for key in category_to_techniques.keys():
			print "=> " +str(key) +" <="
			for tech in category_to_techniques[key].keys():
				print "    -> "+str(tech) +" : " +str(category_to_techniques[key][tech])
				detected += category_to_techniques[key][tech]
		print "=> "+str(detected) +" / " + str(total) +" || [ "+str(float(float(detected)/float(total))*100) +" % ]"

	## return dictionnary
	return category_to_techniques



def generate_techniques_figures(category_to_techniques_to_count):
	##
	## Generate 8 figures from techniques count data
	## category_to_techniques_to_count is a dictionnary obtain by the
	## retrieve_techniques() function.
	##
	## save generated figures in images folder.
	##

	## Define category of analysis
	regression_analysis = ["logistic regression", "cox regression", "multivariate regression",
	"univariate regression", "linear regression", "Poisson regression", "regression tree", "regression analysis"]
	machine_learning_analysis = ["random forest", "artificial neural network", "support vector machine", "machine learning",
	"GLM", "mixed models", "k mean", "clustering", "hierarchical clustering", "decision tree"]
	other_analysis = ["t-test", "Wilcoxon test", "PCA", "chi-squared", "multivariate analysis", "univariate analysis", "discriminant analysis",
	"bioinformatic analysis", "normalization"]

	## Detailed histogramm
	Total_number_of_techniques = 0
	total_regression_techniques = 0
	total_machine_learning_techniques = 0
	total_other_techniques = 0
	global_techniques_data = {}
	for category in category_to_techniques_to_count.keys():

		regression_techniques = 0
		machine_learning_techniques = 0
		other_techniques = 0
		number_of_techniques = 0
		for tech in category_to_techniques_to_count[category].keys():
			number_of_techniques += category_to_techniques_to_count[category][tech]
			Total_number_of_techniques += category_to_techniques_to_count[category][tech]

			if(tech not in global_techniques_data.keys()):
				global_techniques_data[tech] = category_to_techniques_to_count[category][tech]
			else:
				global_techniques_data[tech] += category_to_techniques_to_count[category][tech]

			if(tech in regression_analysis):
				regression_techniques += category_to_techniques_to_count[category][tech]
				total_regression_techniques += category_to_techniques_to_count[category][tech]
			elif(tech in machine_learning_analysis):
				machine_learning_techniques += category_to_techniques_to_count[category][tech]
				total_machine_learning_techniques += category_to_techniques_to_count[category][tech]
			elif(tech in other_analysis):
				other_techniques += category_to_techniques_to_count[category][tech]
				total_other_techniques += category_to_techniques_to_count[category][tech]

		techniques_to_count = category_to_techniques_to_count[category]
		for tech in techniques_to_count.keys():
			techniques_to_count[tech] = float(techniques_to_count[tech]) / float(number_of_techniques) *100


		## Detailed Histogram
		plt.bar(techniques_to_count.keys(), techniques_to_count.values(), color='b', align='center', width=0.3)
		plt.xticks(rotation=45)
		plt.savefig("images/techniques_histogram_"+str(category)+".png")
		plt.close()

		## Pie chart
		regression_techniques = float(regression_techniques) / float(number_of_techniques) * 100
		machine_learning_techniques = float(machine_learning_techniques) / float(number_of_techniques) * 100
		other_techniques = float(other_techniques) / float(number_of_techniques) * 100
		x_vector = [regression_techniques, machine_learning_techniques, other_techniques]
		labels = ["regression techniques", "machine learning techniques", "other techniques"]
		plt.pie(x_vector, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
		plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
		plt.savefig("images/techniques_pie_"+str(category)+".png")
		plt.close()

	## Global histogram
	plt.bar(global_techniques_data.keys(), global_techniques_data.values(), color='b', align='center', width=0.3)
	plt.xticks(rotation=45)
	plt.savefig("images/techniques_histogram_global.png")
	plt.close()

	## Global pie chart
	total_regression_techniques = float(total_regression_techniques) / float(Total_number_of_techniques) * 100
	total_machine_learning_techniques = float(total_machine_learning_techniques) / float(Total_number_of_techniques) * 100
	total_other_techniques = float(total_other_techniques) / float(Total_number_of_techniques) * 100
	x_vector = [total_regression_techniques, total_machine_learning_techniques, total_other_techniques]
	labels = ["total regression techniques", "total machine learning techniques", "total other techniques"]
	plt.pie(x_vector, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
	plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	plt.savefig("images/techniques_pie_global.png")
	plt.close()






"""
TEST SPACE
"""





#get_count_of_articles_type_for_each_disease()

truc = retrieve_techniques(True)
generate_techniques_figures(truc)


"""
print "------[TEST SPACE]------\n"
machin = get_ListOfDirectInteraction("P43403", "P06239")
print "-------------------------------------------------------"
print machin
"""





#plot_pulbications_years("SAVE/run_2/meta")

log_file =  "SAVE/run_14h:1m:25:1/bibotlite.log"
run_folder = "SAVE/run_14h:1m:25:1"
"""
found_cmpt = 0
abstract_list = glob.glob(run_folder+"/abstract/*_abstract.txt")

for abstract in abstract_list:
	data = open(abstract, "r")
	for line in data:
		print line
	data.close()

for abstract_file in abstract_list:
	size = get_cohorte_size(abstract_file)
	if(size != "NA"):
		found_cmpt += 1

print "=> "+str(found_cmpt) +" / "+str(len(abstract_list)) +" [" +str(float(float(found_cmpt)/float(len(abstract_list))*100)) +"]"
"""







"""
diagnostic_abstracts = ["abstract/29363510_abstract.txt", ]
therapeutic_abstracts = ["abstract/29359591_abstract.txt", "25672757_abstract.txt"]
modelisation_abstracts = ["abstract/29333443_abstract.txt"]
abstract = "abstract/25573986_abstract.txt"

import shutil

class_to_count = {"Diagnostic":0, "Therapeutic":0, "Modelisation":0, "Unclassified":0}
for abstract in glob.glob("abstract/*_abstract.txt"):
	theme = get_article_domain(abstract)
	class_to_count[theme] += 1
	pmid = abstract.split("/")
	pmid = pmid[-1]
	shutil.copy(abstract, "articles_classification/"+theme+"/"+pmid)

print class_to_count
"""

#get_all_subjects_from_abstract("abstract")

#describe_articles_type("SAVE/run_2")

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