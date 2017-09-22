"""
Screen the mint database
Molecular interaction database
"""

from Bio import Entrez
from Bio.Entrez import efetch, read
import re
import pprint
import os
from bioservices import *

Entrez.email = 'murlock.raspberypi@gmail.com'


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