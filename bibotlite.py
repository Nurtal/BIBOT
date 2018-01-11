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

	check_date = True
	check_language = True


	##---------------##
	## The Easy Part ##
	##---------------##
	## get meta data on the articles
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


	## Basic check on meta data
	## - check date
	if(int(year) < int(oldest_year_authorized)):
		check_date = False

	## - check language
	if(article_language not in authorized_languages):
		check_language = False


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


	##----------------##
	## The Smart Part ## 
	##----------------##
	## run further analysis on the abstract using nltk
	## => TODO

	## fetch the abstract and convert it to
	## a nltk text object.
	abstract_file_name = "abstract/"+str(pmid)+"_abstract.txt"
	abstract = fetch_abstract(pmid)
	if(abstract):
		save_abstract(abstract, abstract_file_name)
		abstract_text = load_text(abstract_file_name)




##------##
## MAIN ##
##------##




## request
machin = get_ListOfArticles("Big Data AND Sjogren", 15)
machin = get_ListOfArticles("Shen", 150)


for pmid in machin:

	evaluate_article(pmid)

	"""

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