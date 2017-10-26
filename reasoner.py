
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


def create_text():
	## Create nltk Text from text file
	## write the text file from abstract
	## fetched by BIBOT.
	## Return a nltk text object

	## -> Create the input file text
	## from pubmed abstract file.
	## -> Get text
	abstract_file = open("fetched/pubmed_abstract.txt", "r")
	text = ""
	for line in abstract_file:
		line = line.replace("\n", "")
		if(line[0] != ">"):
			text += line
	abstract_file.close()

	# -> write input file
	text_file = open("fetched/nltk_text.txt", "w")
	text_file.write(text)
	text_file.close()

	## -> Create the Text object from the input
	## text file
	text_file=open('fetched/nltk_text.txt','rU')
	raw=text_file.read()
	tokens = nltk.word_tokenize(raw)
	text = nltk.Text(tokens)

	## Return the nltk text object
	return text




### TEST SPACE ###
stuff = create_text()
print stuff
text1.concordance("monstrous")