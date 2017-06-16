"""
Crawler
"""

import rules_manager
import bibliosearch


def gather_doc_rules():
	"""
	[IN PROGRESS]
	"""


	rules = rules_manager.get_all_rules_from_file("HLA_phase1_all_rules_filtered.txt")
	cmpt = 0
	for rule in rules.values():

		result_file = open("results/abstract_for_rule_"+str(cmpt)+".txt", "w")
		query = rules_manager.generate_querry_from_rule(rule)
		pmid_list = bibliosearch.get_ListOfArticles(query, 1500)
		for pmid in pmid_list:
			try:
				test = bibliosearch.fetch_abstract(pmid)
				print "[+] Read "+str(pmid)
				print test
				result_file.write(">PMID;"+str(pmid)+"\n")
				result_file.write(test+"\n")
			except:
				result_file.write("[!] Can't read "+str(pmid)+"\n")
				print "[!] Can't read "+str(pmid)
		result_file.close()


### TEST SPACE ###
gather_doc_rules()