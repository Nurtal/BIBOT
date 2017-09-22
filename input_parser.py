"""
Set of functions to parse
input from the user.
"""

import ncbi_crawler
import mygene

def search_parser(query, query_type):
	"""
	-> In progress
	"""

	##-----------------------------------##
	## Deal with a classic biblio search ##
	##-----------------------------------##
	if(query_type == 1):
		
		## get list of PMID
		list_of_articles = ncbi_crawler.get_ListOfArticles(query, 10)

		## Write results in a data file
		output_file = open("fetched/pubmed_abstract.txt", "w")
		for pmid in list_of_articles:
			abstract = ncbi_crawler.fetch_abstract(pmid)
			abstract = abstract.encode('ascii', 'ignore')
			output_file.write(">"+str(pmid)+"\n")
			output_file.write(abstract+"\n")
		output_file.close()


	##-----------------------##
	## Get Gene informations ##
	##-----------------------##
	elif(query_type == 2):

		## Collect information
		mg = mygene.MyGeneInfo()
		gene_info = mg.query(query, size=1)
		gene_name = gene_info["hits"][0]["name"]
		gene_symbol = gene_info["hits"][0]["symbol"]
		gene_taxid = gene_info["hits"][0]["taxid"]
		gene_entrez = gene_info["hits"][0]["entrezgene"]
		gene_id = gene_info["hits"][0]["_id"]

		## Write results in a data file
		## <choucroute> separateur
		output_file = open("fetched/gene_information.csv", "w")
		output_file.write("name<choucroute>"+str(gene_name)+"\n")
		output_file.write("symbol<choucroute>"+str(gene_symbol)+"\n")
		output_file.write("taxid<choucroute>"+str(gene_taxid)+"\n")
		output_file.write("entrezgene<choucroute>"+str(gene_entrez)+"\n")
		output_file.write("id<choucroute>"+str(gene_id)+"\n")
		output_file.close()

	##------------------##
	## Work in progress ##
	##------------------##
	elif(query_type == 3):
		print "Work in progress"

	##---------------------------------##
	## Wrong parameters initialization ##
	##---------------------------------##
	else:
		print "What do you looking for ?"
