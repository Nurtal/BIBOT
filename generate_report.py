"""
Stuff to generate report
"""




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