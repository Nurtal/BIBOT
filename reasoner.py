"""
testing stuff
"""


from lxml import etree

def get_number_of_relations(file):
	"""
	Return the number (int) of relations describe
	in file
	"""
	tree = etree.parse("template.xml")
	for number in tree.xpath("/object/number_of_relations"):
		num = int(number.text)
	return num


def get_relations(file):
	"""
	Return the list of relations
	descrbibe in file
	"""
	number_of_relations = get_number_of_relations(file)

	list_of_relations = []
	tree = etree.parse(file)
	relation_number = 0
	
	for relation in xrange(number_of_relations):
		relation_number = relation_number + 1
		relation_in_str = ""
		left_members = ""
		operation = ""
		right_members = ""
		left_count = 0
		right_count = 0
		for item_left in tree.xpath("/object/relation_"+str(relation_number)+"/item_left"):
			left_count = left_count + 1
			if left_count > 1:
				left_members = left_members + " AND " + item_left.text
			else:
				left_members = left_members + item_left.text
		for action in tree.xpath("/object/relation_"+str(relation_number)+"/action"):
			operation = operation + action.text
		for item_right in tree.xpath("/object/relation_"+str(relation_number)+"/item_right"):
			right_count = right_count + 1
			if right_count > 1:
				right_members = right_members + " AND " + item_right.text
			else:
				right_members = right_members + item_right.text
		relation_in_str = left_members + " " +operation+ " " + right_members
		list_of_relations.append(relation_in_str)

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

def get_parent(file):
	"""
	Return the list of parent
	descrbibe in file
	"""
	list_of_parent = []
	tree = etree.parse(file)
	for parent in tree.xpath("/object/parent"):
		list_of_parent.append(parent.text)
	return list_of_parent



def get_all_relations(file, list_of_relations):
	"""
	get all relations, including parent's relation
	(recusrive function)

	to test with the new syntax
	"""
	relation_list = get_relations(file)
	list_of_relations = list_of_relations + relation_list
	parent_list = get_parent(file)
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

test_sentence_1 = "<item_1>monocyte</item1> <action>increase</action> <item_2>lymphocyte</item_2>"
test_sentence_2 = "<item_1>monocyte</item1> <action>increase</action> <item_2>lymphocyte</item_2>"

relations = get_relations("monocyte.xml")
print relations

nature = get_nature("monocyte.xml")
print nature

parent = get_parent("monocyte.xml")
print parent

print "--------------------------"

machin = get_all_relations("monocyte.xml", relations)
print machin
