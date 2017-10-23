"""
Stuff to generate report
"""

import os.path
import webbrowser



def write_abstract_report():
	"""
	IN PROGRESS
	-> Classic report generation,
	write an html file with the abstract fetched from PUBMED

	-> TODO:
		- write title (query)
		- Write title of the article from PMID
		- Create index table for articles
	"""

	report = open("reports/biblio.html", "w")
	report.write("<html>\n")
	report.write("<title>BIBOT - biblio report</title>\n")

	## Access Data
	data_file = open("fetched/pubmed_abstract.txt", "r")
	for line in data_file:
		line = line.replace("\n", "")
		if(line[0] == ">"):
			report.write("<h2>"+line+"</h2>\n")
		else:
			report.write("<p>"+line+"</p>\n")
	data_file.close()
	report.write("</html>\n")

	report.close()


def write_abstract_report_dev():
	"""
	IN PROGRESS
	-> Classic report generation,
	write an html file with the abstract fetched from PUBMED

	-> TODO:
		- write title (query)
		- Write title of the article from PMID
		- Create index table for articles
	"""

	report = open("reports/biblio.html", "w")
	report.write("<html>\n")
	report.write("<body background=\"report_background.jpg\">\n")
	report.write("<title>BIBOT - biblio report</title>\n")


	## Gene Section
	if(os.path.isfile("fetched/gene_information.csv")):
		
		## Access Data
		gene_values_of = {}
		gene_information = open("fetched/gene_information.csv", "r")
		for line in gene_information:
			line = line.replace("\n", "")
			line_in_array = line.split("<choucroute>") # Yep
			gene_values_of[str(line_in_array[0])] = str(line_in_array[1])
		gene_information.close()

		## Write Section
		report.write("<h1>Gene Information</h1>\n")

		report.write("<table style=\"width:100%\">\n")
		
		report.write("<tr>\n")
		report.write("<td>name</td>\n")
		report.write("<td>"+gene_values_of["name"]+"</td>\n")
		report.write("</tr>\n")
		
		report.write("<tr>\n")
		report.write("<td>Symbol</td>\n")
		report.write("<td>"+gene_values_of["symbol"]+"</td>\n")
		report.write("</tr>\n")

		report.write("<tr>\n")
		report.write("<td>ID</td>\n")
		report.write("<td>"+gene_values_of["id"]+"</td>\n")
		report.write("</tr>\n")

		report.write("<tr>\n")
		report.write("<td>Entrez Gene ID</td>\n")
		report.write("<td>"+gene_values_of["entrezgene"]+"</td>\n")
		report.write("</tr>\n")

		report.write("<tr>\n")
		report.write("<td>Taxid</td>\n")
		report.write("<td>"+gene_values_of["taxid"]+"</td>\n")
		report.write("</tr>\n")

		report.write("</table>\n")

	## Protein Section
	if(os.path.isfile("fetched/protein_information.csv")):
		
		## Access Data
		prot_values_of = {}
		prot_information = open("fetched/protein_information.csv", "r")
		for line in prot_information:
			line = line.replace("\n", "")
			line_in_array = line.split(",")
			prot_values_of[str(line_in_array[0])] = str(line_in_array[1])
		prot_information.close()

		## Write Section
		report.write("<h1>Protein Information</h1>\n")

		report.write("<table style=\"width:100%\">\n")
		
		report.write("<tr>\n")
		report.write("<td>UniProt ID</td>\n")
		report.write("<td>"+prot_values_of["uniprot_id"]+"</td>\n")
		report.write("</tr>\n")

		report.write("</table>\n")


	## Access Data
	data_file = open("fetched/pubmed_abstract.txt", "r")
	for line in data_file:
		line = line.replace("\n", "")
		if(line[0] == ">"):
			report.write("<h2>"+line+"</h2>\n")
		else:
			report.write("<p>"+line+"</p>\n")
	data_file.close()
	report.write("</html>\n")

	report.close()

	## Display report
	##webbrowser.open("reports/biblio.html") # not working
