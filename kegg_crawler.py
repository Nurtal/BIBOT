"""
Screen the KEGG database
"""

from Bio import Entrez
from Bio.Entrez import efetch, read
import re
import pprint
import os
from bioservices import *

Entrez.email = 'murlock.raspberypi@gmail.com'


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


def get_involved_pathways(gene_id):
	"""
	-> Create a connection to the KEGG database,
	   return a list of pathways (with id) where the
	   gene_id (entrez stuff) is involved
	-> "hsa" is the human reference for the KEGG database
	"""
	k = KEGG(verbose=False)
	pathways = k.get_pathway_by_gene(str(gene_id), "hsa")

	return pathways


def show_pathway(pathway_name):
	"""
	-> Display pathway using KEGG database
	-> Retrieve gene id in gene_information.csv
	-> Retrieve pathway id from pathway name in pathways_involved.csv
	IN PROGRESS
	"""
	k = KEGG(verbose=False)

	pathway_id = "undef"
	gene_id = "undef"

	## Get pathway ID
	pathway_file = open("fetched/pathways_involved.csv", "r")
	for line in pathway_file:
		line = line.replace("\n", "")
		line_in_array = line.split(",")
		if(line_in_array[1] == pathway_name):
			pathway_id = line_in_array[0]
	pathway_file.close()

	## Get gene ID
	gene_file = open("fetched/gene_information.csv", "r")
	for line in gene_file:
		line = line.replace("\n", "")
		line_in_array = line.split("<choucroute>") # Yep
		if(line_in_array[0] == "entrezgene"):
			gene_id = line_in_array[1]
	gene_file.close()

	## Show the figure
	if(gene_id != "undef" and pathway_id != "undef"):
		k.show_pathway(str(pathway_id), keggid={str(gene_id): "red"})


	




##---------------------------------------------------##
## Exemple d'utilisation des fonctions de l'API KEGG ##
##---------------------------------------------------##
#k = KEGG(verbose=False)
#k.find("hsa", "zap70")
#pathway = k.get_pathway_by_gene("7535", "hsa")
#print pathway
#k.show_pathway("hsa04064", keggid={"7535": "red"})