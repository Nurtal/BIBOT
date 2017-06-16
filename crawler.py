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
	for rule in rules.values():
		query = generate_querry_from_rule(rule)
		pmid_list = bibliosearch.get_ListOfArticles(query, 1500)
		for pmid in machin:
			try:
				test = fetch_abstract(pmid)
				print "["+str(pmid)+"]\n"
				print test
			except:
				print "[!] Can't read "+str(pmid)

### TEST SPACE ###
gather_doc_rules()