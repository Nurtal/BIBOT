#!/usr/bin/python
#coding: utf8 
from __future__ import unicode_literals
##==================================================================##
##=======================>BIBOTLIGHT<===============================##
##==================================================================##
## -> Light version of BIBOT, focus on pubmed fetching and natural  ## 
## langage processing for the selected articles.                    ##
##
##
##

##-------------##
## IMPORTATION ##
##-------------##
from Bio import Entrez
from Bio.Entrez import efetch, read
from unidecode import unidecode
import nltk
import itertools
import os



##------------##
## PARAMETERS ##
##------------##
Entrez.email = 'murlock.raspberypi@gmail.com'




##-----------##
## FUNCTIONS ##
##-----------##


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
	=> Not working
	"""
	handle = efetch(db='pubmed', id=pmid, retmode='xml', )
	informations = read(handle)

	stuff = informations[u'PubmedArticle'][0] 
	date = stuff[u'PubmedData']["History"][1]
	month = date[u'Month']
	day = date[u'Day']
	year = date[u'Year']

	print month
	print day
	print year

	return "choucroute"


def get_ListOfArticles(term, max_count):
	"""
	return the list of pmid article conatining
	the term.
	"""
	h = Entrez.esearch(db='pubmed', retmax=max_count, term=term)
	result = Entrez.read(h)
	listOfArticles = result["IdList"]

	return listOfArticles;



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


def load_text(text_file):
	##
	## -> Create  and return nltk Text 
	## object from text_file
	##

	## -> Create the Text object from the input
	## text file
	text_file=open(text_file,'rU')
	raw=text_file.read()
	raw = raw.decode('utf8')
	tokens = nltk.word_tokenize(raw)
	text = nltk.Text(tokens)

	## Return the nltk text object
	return text


def evaluate_article(pmid):
	##
	## [IN PROGRESS]
	##
	## -> Test if the abstract is cool
	## -> return true or false
	##

	##------------------------##
	## Parameters for filters ##
	##------------------------##
	oldest_year_authorized = 2008
	authorized_languages = [u'eng']

	valid_article = False
	check_date = True
	check_language = True
	validation_check_keywords_1 = False
	validation_check_keywords_2 = False



	##---------------##
	## The Easy Part ##
	##---------------##
	## get meta data on the articles
	try:
		handle = efetch(db='pubmed', id=pmid, retmode='xml', )
		informations = read(handle)
		stuff = informations[u'PubmedArticle'][0] 
		
		## get date from the history attribute, select
		## the date of acceptation.
		date = stuff[u'PubmedData']["History"][1]
		month = date[u'Month']
		day = date[u'Day']
		year = date[u'Year']

		## get the name of the review
		journal_name = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'MedlineJournalInfo'][u'MedlineTA']
		
		## get the keywords for the articles
		## the format is a bit strange, may have to be carreful
		## with this data (mix of strings and unicode elements)
		keywords_list = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'KeywordList']

		## Get the author's conflict of interest,
		## because we can.
		try:
			conflict_of_interest = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'CoiStatement']
		except:
			conflict_of_interest = "NA"

		## Get title of the article
		article_title = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'Article'][u'ArticleTitle']

		## Get language of the article
		article_language = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'Article'][u'Language'][0]

	except:
		return (False,False,False)

	##----------------##
	## The Smart Part ## 
	##----------------##
	## run further analysis on the abstract using nltk

	## fetch the abstract and convert it to
	## a nltk text object.
	abstract_file_name = "abstract/"+str(pmid)+"_abstract.txt"
	abstract = fetch_abstract(pmid)
	if(abstract):
		save_abstract(abstract, abstract_file_name)
		abstract_text = load_text(abstract_file_name)
		
		## Play with tokenization and chunking
		## Get all the commun names in the abstract
		names_found_in_abstract = []
		try:
			tokens = nltk.word_tokenize(abstract.encode('utf8'))
			tagged = nltk.pos_tag(tokens)
			entities = nltk.chunk.ne_chunk(tagged)
		except:
			print "[WARNINGS] => can't perform nlp operation"
			entities = []

		for item in entities:
			try:
				if(item[1] in ["NN", "NNS", "NNP"]):
					if(item[0] not in names_found_in_abstract):
						names_found_in_abstract.append(item[0])
			except:
				## Somethig went wrong
				choucroute = True
				
		## -> Biology keywords check
		## -> Artificial intelligence keywords check
		IA_keywords = ["algorithm", "machine" "learning", "neural", "network", "statistic", "deep"]
		Clinical_keywords = ["Sjogren" "lupus", "autoimmunity"]
		for item in names_found_in_abstract:
			if(item in IA_keywords):
				validation_check_keywords_1 = True
			if(item in Clinical_keywords):
				validation_check_keywords_2 = True
		
	##--------------##
	## PASS OR FAIL ##
	##--------------##
	## General check phase
	easy_check_passed = False
	smart_check_passed = False

	## Basic check on meta data
	## - check date
	if(int(year) < int(oldest_year_authorized)):
		check_date = False

	## - check language
	if(article_language not in authorized_languages):
		check_language = False

	## Easy Filter
	if(check_date and check_language):
		easy_check_passed = True

	## Complex filter
	if(validation_check_keywords_1 and validation_check_keywords_2):
		smart_check_passed = True

	## Global check
	if(easy_check_passed and smart_check_passed):
		valid_article = True

	##-------------##
	## SAVING DATA ##
	##-------------##
	## Write and delete files
	if(valid_article):

		## Save meta data in a text file
		## for further use
		title_line = u'>Title;'+unicode(article_title)+u"\n"
		date_line = u'>Date;'+unicode(day)+u"/"+unicode(month)+u"/"+unicode(year)+u"\n"
		journal_line = u">Journal;"+unicode(journal_name)+u"\n"
		conflict_of_interest_line = u">Conflict;"+unicode(conflict_of_interest)+u"\n"
		meta_data = open("meta/"+str(pmid)+".csv", "w")
		meta_data.write(title_line.encode('utf8'))
		meta_data.write(date_line.encode('utf8'))
		meta_data.write(journal_line.encode('utf8'))
		meta_data.write(conflict_of_interest_line.encode('utf8'))
		meta_data.close()

	else:
		## Delete the abstract
		try:
			if(abstract):
				os.remove(abstract_file_name)
		except:
			print "[WARNING] => can't delete "+str(abstract_file_name)

	##------------------##
	## RETURN SOMETHING ##
	##------------------##
	## return True if the article pass the 
	## evaluation, else False.
	return (valid_article, easy_check_passed, smart_check_passed)



def get_huge_list_of_artciles(keywords):
	##
	## Create all possible comination of at least two elements
	## from the keywords list. then use these combination
	## to create request (using only the AND operator for now)
	## and screen the pubmed database.
	##
	## return the list of all articles found
	##

	## init variables
	huge_list_of_PMID = []

	## make all combination of at least 2 item in keywords
	combination_list = []
	for x in xrange(2, len(keywords)):
		machin = itertools.combinations(keywords, x)
		for truc in machin:
			combination_list.append(list(truc))

	## create request
	for items_set in combination_list:
		request = ""
		for item in items_set:
			request += item +" AND "
		request = request[:-5]
		
		## screening pubmed
		results_PMID = get_ListOfArticles(request, 9999999)
		
		## increment huge list of pmid
		for pmid in results_PMID:
			if(pmid not in huge_list_of_PMID):
				huge_list_of_PMID.append(pmid)

	return huge_list_of_PMID



##------##
## MAIN ##
##------##




## request
#machin = get_ListOfArticles("Big Data AND Sjogren", 1)
#evaluate_article(machin[0])


request_term = ["big data", "artificial intelligence", "autoimmunity", "Sjogren", "RA", "SLE", "lupus"]
truc = get_huge_list_of_artciles(request_term)
total_number = len(truc)
fetched = 0
first_fiter_passed = 0
last_filter_passed = 0
cmpt = 0
for article in truc:
	valid = evaluate_article(article)
	filter_1_status = "FAILED"
	filter_2_status = "FAILED"
	if(valid[0]):
		fetched += 1
	if(valid[1]):
		first_fiter_passed += 1
		filter_1_status = "PASSED"
	if(valid[2]):
		last_filter_passed += 1
		filter_2_status = "PASSED"
	cmpt += 1

	#print str(cmpt) +" || " +str(fetched) + " || " +str(total_number)
	print "|| "+str(cmpt) +" [PROCESSED] || "+ str(fetched) + " [SELECTED] || FIRST FILTERS ["+filter_1_status+ "] || LAST FILTER ["+filter_2_status+ "] || "+str(float((float(cmpt)/float(total_number))*100)) + "% [COMPLETE]"


print "FIRST FILTER PASS => "+str(first_fiter_passed)
print "LAST FILTER PASS => "+str(last_filter_passed)
"""
## tokenization exemple
from nltk.corpus import treebank
sentence = "the dog is green and the rabbit is mean"
tokens = nltk.word_tokenize(sentence)
tagged = nltk.pos_tag(tokens)
entities = nltk.chunk.ne_chunk(tagged)
print entities
entities.draw()
"""



"""
for pmid in machin:
	evaluate_article(pmid)

	## get abstract
	abstract_file_name = str(truc)+"_abstract.txt"
	abstract = fetch_abstract(truc)
	save_abstract(abstract, abstract_file_name)
	text = load_text(abstract_file_name)

	## nltk playground
	## with abstract
	print len(text)
	text.concordance("data")
	fdist1 = nltk.FreqDist(text)
	print fdist1
	print fdist1.most_common(20)
	print text.collocations()


	"""