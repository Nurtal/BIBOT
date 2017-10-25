
##
## Play with natural langage
##


import nltk
from nltk.corpus.reader.plaintext import PlaintextCorpusReader
#nltk.download()
#from nltk.book import *


def create_corpus():
	## Create corpus from abstract
	## fetched by BIBOT
	## return a corpus object

	## Read the abstract result file
	abstract_to_content = {}
	abstract_file = open("fetched/pubmed_abstract.txt", "r")
	for line in abstract_file:
		line = line.replace("\n", "")
		if(line[0] == ">"):
			abstract = line[1:]
			abstract_to_content[abstract] = ""
		else:
			content = line
			abstract_to_content[abstract] = content
	abstract_file.close()

	## create files
	for key in abstract_to_content.keys():
		text_file = open("fetched/corpus/"+str(key)+".txt", "w")
		text_file.write(abstract_to_content[key])
		text_file.close()

	## ntlk magical lines
	corpusdir = 'fetched/corpus/'
	newcorpus = PlaintextCorpusReader(corpusdir, '.*')

	return newcorpus


truc = create_corpus()
print truc.raw().strip()