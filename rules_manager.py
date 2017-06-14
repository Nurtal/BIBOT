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



### TEST SPACE ###
machin = get_all_rules_from_file("HLA_phase1_all_rules_filtered.txt")
print machin