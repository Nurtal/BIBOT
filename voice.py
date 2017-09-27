"""
Text to speach module
using pyttsx
work on python 2.7 32 bit (under windows, not tested under linux yet)
"""

import pyttsx
import random


english_voice = "HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Speech\Voices\Tokens\TTS_MS_EN-US_ZIRA_11.0"
french_voice = "HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Speech\Voices\Tokens\TTS_MS_FR-FR_HORTENSE_11.0"

def say_something():
	"""
	test function
	"""
	engine = pyttsx.init()
	engine.say('Franchement, Precise Sad se ne vaut pas une bonne choucroute')
	engine.runAndWait()




def say_biological_truth():
	"""
	A little prank
	"""

	## Collection of sentences
	engine = pyttsx.init()
	biological_truth = [
		"Le probleme avec la biologie, c'est les biologistes",
		"Je veux faire un test statistique avec cinquante et une variable pour trois patients.",
		"Nous avons quatre vingt dix pourcent de valeurs manquantes. Nous allons les infairai."
	]

	## Pick and say
	selection = random.randint(0,len(biological_truth)-1)
	engine.say(biological_truth[selection])
	engine.runAndWait()



def write_speech():
	"""
	-> Get general informations about the request
	-> Write the text to say
	[IN PROGRESS]
	"""

	## Check if a search has been done


	## General stuff about the search
	## looking for a gene ? abstract ?

	## Count the number of pathway involved




## TEST SPACE ##

say_biological_truth()
#say_something()
