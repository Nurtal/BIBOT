"""
Functions for parsing generated
rules, extract information and so on...
Prepare query for the web crawling part
"""

def get_all_rules_from_file(data_file_name):
	"""
	-> Store all rules present in data_file_name
	into a dict structure
	-> return the dict structure
	   with:
	   	-left terms (list)
	   	-right terms (list)
	   	-support (float, cast in string)
	   	-confidence (float, cast in string)
	   	-lift (float, cast in string)
	"""

	rules_structure = {}
	input_data = open(data_file_name, "r")
	cmpt = 0
	for line in input_data:
		line = line.split("\n")
		line = line[0]

		if(cmpt > 0):
			line_in_array = line.split("=>")
			left_part = line_in_array[0]
			right_part = line_in_array[1]

			## Isolate the left part of the equation
			left_part = left_part.split("{")
			left_part = left_part[1]
			left_part = left_part.replace("}", "")
			left_terms = left_part.split(",")

			## Isolate the right part of the equation
			right_part = right_part.split("\";")
			stat_part = right_part[1]
			right_part = right_part[0]
			right_part = right_part.replace("}", "")
			right_part = right_part.replace("{", "")
			right_terms = right_part.split(",")

			stat_part = stat_part.split(";")
			support = stat_part[0]
			confidence = stat_part[1]
			lift = stat_part[2]

			## Fill the dictionnary
			rules_structure[cmpt] = {"left_terms":left_terms,
									 "right_terms":right_terms,
									 "support":support,
									 "confidence":confidence,
									 "lift":lift}

		cmpt+=1
	input_data.close()
	return rules_structure




def generate_querry_from_rule(rule):
	"""
	-> Use the left and right members of a rule
	   to generate a query in order to fetch abstracts
	-> rule is a dict structure generate by the get_all_rules_from_file()
	   function
	-> return the query (a string)
	"""

	## Structure initialisation
	query_left_part = ""
	query_left_elements = []
	query_right_part = ""
	query_right_elements = []
	final_query = ""

	## Get information
	left_terms = rule['left_terms']
	right_terms = rule['right_terms']
	
	## Create first part of the query, from
	## left termes of the rule
	for element in left_terms:
		element_in_array = element.split("=")
		terms = element_in_array[0]
		terms_in_array = terms.split(".")
		terms_in_array = terms_in_array[1:] # drop the "X"

		for term in terms_in_array:
			if(term not in query_left_elements):
				query_left_elements.append(term)


	## Write the first part of the query
	for element in query_left_elements:
		query_left_part += str(element) +" AND "
	query_left_part = query_left_part[:-5]


	## Create second part of the query, from
	## right terms of the rule
	for element in right_terms:
		element_in_array = element.split("=")
		terms = element_in_array[0]

		terms_in_array = terms.split(".")
		terms_in_array = terms_in_array[1:] # drop the "X"

		for term in terms_in_array:
			if(term not in query_right_elements and term not in query_left_elements):
				query_right_elements.append(term)


	## Write the first part of the query
	for element in query_right_elements:
		query_right_part += str(element) +" OR "
	query_right_part = query_right_part[:-4]

	## Merge the 2 parts of the query
	if(len(query_right_part) > 0):
		final_query = query_left_part +" OR "+query_right_part

	return final_query




### TEST SPACE ###
#machin = get_all_rules_from_file("HLA_phase1_all_rules_filtered.txt")
#test_rule = machin[5]
#generate_querry_from_rule(test_rule)
