"""
testing stuff
"""


from lxml import etree




def get_relations(file):
	"""
	Return the list of relations
	descrbibe in file
	"""
	list_of_relations = []
	tree = etree.parse(file)
	for relation in tree.xpath("/object/relation"):
		list_of_relations.append(relation.text)
	return list_of_relations


def get_nature(file):
	"""
	Return the list of nature
	descrbibe in file
	"""
	list_of_nature = []
	tree = etree.parse(file)
	for nature in tree.xpath("/object/nature"):
		list_of_nature.append(nature.text)
	return list_of_nature


def get_all_relations(file, list_of_relations):
	"""
	get all relations, including parent's relation
	(recusrive function)
	"""
	relation_list = get_relations(file)
	list_of_relations = list_of_relations + relation_list
	parent_list = get_nature(file)
	for parent in parent_list:
		if(parent != "primal"):
			newFile = parent +".xml"
			list_of_relations = get_all_relations(newFile, list_of_relations)
	
	return list_of_relations


"""
TEST SPACE
"""

relationshipt_dict = {"increase":[{"equivalent":["augmentation", "developpement"]}, {"opposite":["decrease", "decline"]}]}

sentence = "monocyte increase lymphocyte"


list_of_relations = []
list_of_relations = get_all_relations("monocyte.xml", list_of_relations)
for elt in list_of_relations:
	print elt + " & " + sentence