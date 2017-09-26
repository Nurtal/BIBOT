"""
Set of functions to parse
input from the user.
"""

import ncbi_crawler
import kegg_crawler
import mygene
import settings

def search_parser(query, query_type):
	"""
	[IN PROGRESS]
	"""

	##-----------------------------------##
	## Deal with a classic biblio search ##
	##-----------------------------------##
	if(query_type == 1):
		
		## Get the parameters for the search
		parameters = settings.read_settings()

		## get list of PMID
		list_of_articles = ncbi_crawler.get_ListOfArticles(query, parameters['max_abstract_return'])

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

		## Check if we found something
		if(len(gene_info["hits"]) > 0):
			gene_name = gene_info["hits"][0]["name"]
			gene_symbol = gene_info["hits"][0]["symbol"]
			gene_taxid = gene_info["hits"][0]["taxid"]
			gene_entrez = gene_info["hits"][0]["entrezgene"]
			gene_id = gene_info["hits"][0]["_id"]

			##-------------------------##
			## Get Protein Information ## 
			##-------------------------##
			write_protein_information(gene_entrez)

			##-----------------------##
			## Get Pathways Involved ##
			##-----------------------##
			pathways_fetched = kegg_crawler.get_involved_pathways(gene_entrez)
			pathways_file = open("fetched/pathways_involved.csv", "w")
			for pathway in pathways_fetched:
				pathways_file.write(pathway+","+pathways_fetched[pathway]+"\n")
			pathways_file.close()


		else:
			gene_name = "No Gene Found"
			gene_symbol = "Nothing Found"
			gene_taxid = "Nothing Found"
			gene_entrez = "Nothing Found"
			gene_id = "Nothing Found"

		## Write results in a data file
		## <choucroute> separateur, beacause random
		## char always present in gene name
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



def write_protein_information(entrez_gene_id):
	"""
	IN PROGRESS
	"""
	from PyEntrezId import Conversion

	# include your email address
	Id = Conversion('dummyemail@dummybunny.info')
	UniProtId = Id.convert_entrez_to_uniprot(entrez_gene_id)
	
	## Write results in file
	output_file = open("fetched/protein_information.csv", "w")
	output_file.write("uniprot_id,"+str(UniProtId)+"\n")
	output_file.close()


def show_pathway_involved(pathway_name):
	"""
	-> Wrapper around the show_pathway function
	   from the kegg_crawler file
	"""
	kegg_crawler.show_pathway(pathway_name)
