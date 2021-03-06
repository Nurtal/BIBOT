"""
Graphical Interface for BIBOT
"""


from Tkinter import *
#import ImageTk
from PIL import ImageTk, Image
from bibliosearch import *
import input_parser
import generate_report
import settings
import webbrowser
import voice
import os.path
import os


"""Test Stuff"""

def testFunction(machin):
	print "Testing "+str(machin)


def open_uniprot_page():
	"""
	-> open the uniprot entry (web page)
	   for the uniprot_id
	"""

	## Get the uniprot ID
	uniprot_id = "undef"
	protein_data_file = open("fetched/protein_information.csv", "r")
	for line in protein_data_file:
		line = line.replace("\n", "")
		line_in_array = line.split(",")
		if(line_in_array[0] == "uniprot_id"):
			uniprot_id = line_in_array[1]
	protein_data_file.close()

	## open page
	webbrowser.open("http://www.uniprot.org/uniprot/"+str(uniprot_id))


def open_mint_page():
	"""
	-> open the mint entry (web page)
	   for the uniprot_id
	"""

	## Get the uniprot ID
	uniprot_id = "undef"
	protein_data_file = open("fetched/protein_information.csv", "r")
	for line in protein_data_file:
		line = line.replace("\n", "")
		line_in_array = line.split(",")
		if(line_in_array[0] == "uniprot_id"):
			uniprot_id = line_in_array[1]
	protein_data_file.close()

	## open page
	url_address = "http://mint.bio.uniroma2.it/cgi-bin/advResults.py?proteinGeneName=0&protein="+str(uniprot_id)+"&taxon=Homo%20sapiens&interaction=0&interactionType=0&methods=0&uniprot=0&pmid=0http://mint.bio.uniroma2.it/cgi-bin/advResults.py?proteinGeneName=0&protein=P01111&taxon=Homo%20sapiens&interaction=0&interactionType=0&methods=0&uniprot=0&pmid=0"
	webbrowser.open(url_address)



def search(query, query_type):
	"""
	Wrapper around the search parser
	from input_parser
	"""
	input_parser.search_parser(query, query_type)
	button_report_general['state'] = "normal"

	## Get Biblio information
	if(query_type == 1):

		## get number of abstract found
		number_of_abstract_returned = 0
		abstract_file_fetched = open("fetched/pubmed_abstract.txt", "r")
		for line in abstract_file_fetched:
			line = line.replace("\n", "")
			if(line[0] == ">"):
				number_of_abstract_returned += 1
		abstract_file_fetched.close()
		label_abstract_found["text"] = "Abstract Found : "+str(number_of_abstract_returned)

	## Get gene information
	if(query_type == 2):

		data_file = open("fetched/gene_information.csv", "r")
		for line in data_file:
			line = line.replace("\n", "")
			line_in_array = line.split("<choucroute>") # Never mind

			## Update gene label
			if(line_in_array[0] == "name"):
				label_gene_name["text"] = "Name : "+str(line_in_array[1])
			elif(line_in_array[0] == "symbol"):
				label_gene_symbol["text"] = "Symbol : "+str(line_in_array[1])
			elif(line_in_array[0] == "taxid"):
				label_gene_taxid["text"] = "Taxid : "+str(line_in_array[1])
			elif(line_in_array[0] == "entrezgene"):
				label_gene_entrez["text"] = "Entrez-Gene : "+str(line_in_array[1])
			elif(line_in_array[0] == "id"):
				label_gene_id["text"] = "id : "+str(line_in_array[1])
		data_file.close()

		## Clean list of current Pathways involved
		listbox_pathways.delete(0, END)

		## Add Pathways involved information
		pathways_file = open("fetched/pathways_involved.csv", "r")
		for line in pathways_file:
			line = line.replace("\n", "")
			line_in_array = line.split(",")

			## Update pathways listbox
			listbox_pathways.insert(END, line_in_array[1])

		pathways_file.close()

		## Add protein information
		protein_data_file = open("fetched/protein_information.csv", "r")
		for line in protein_data_file:
			line = line.replace("\n", "")
			line_in_array = line.split(",")

			## Update protein label
			## Update buttun state
			if(line_in_array[0] == "uniprot_id"):
				label_uniprotId["text"] = "UniProt ID : "+str(line_in_array[1])
				button_uniprot["state"] = "normal"
				button_mint["state"] = "normal"
		protein_data_file.close()


def display_pathway(pathway_name):
	"""
	-> Wrapper around the show_pathway_involved function
	  from the input_parser file
	"""
	input_parser.show_pathway_involved(pathway_name)



def reset_search():
	"""
	IN PROGRESS
	"""

	## Reset Gene Information
	label_gene_name["text"] = "Name : NA"
	label_gene_symbol["text"] = "Symbol : NA"
	label_gene_taxid["text"] = "Taxid : NA"
	label_gene_entrez["text"] = "Entrez-Gene : NA"
	label_gene_id["text"] = "id : NA"

	## Delete gene info file
	if(os.path.exists("fetched/gene_information.csv")):
		os.remove("fetched/gene_information.csv")

	## Clean list of current Pathways involved
	listbox_pathways.delete(0, END)

	## Delete Pathway file
	if(os.path.exists("fetched/pathways_involved.csv")):
		os.remove("fetched/pathways_involved.csv")

	## Reset Protein information
	label_uniprotId["text"] = "UniProt ID : NA"
	button_uniprot["state"] = "disabled"
	button_mint["state"] = "disabled"

	## Delete Protein file
	if(os.path.exists("fetched/protein_information.csv")):
		os.remove("fetched/protein_information.csv")

	## Reset abstract info
	label_abstract_found["text"] = "Abstract Found : NA"

	## Delete abstract file
	if(os.path.exists("fetched/pubmed_abstract.txt")):
		os.remove("fetched/pubmed_abstract.txt")





"""Interface"""

root = Tk()
#img = ImageTk.PhotoImage(Image.open("bibot_logo.png"))
root.title('BIBOT Interface')
root['bg']='white'


## main frame 

# frame 1
Frame1 = Frame(root, borderwidth=2, relief=GROOVE)
Frame1.grid(column=0, row=1, padx=25, pady=25)

# frame 2
Frame2 = Frame(root, borderwidth=2, relief=GROOVE)
Frame2.grid(column=1, row=1, padx=25, pady=25)

# frame 3
Frame3 = Frame(root, borderwidth=2, relief=GROOVE)
Frame3.grid(column=2, row=1, padx=25, pady=25)

## Frame 4
Frame4 = Frame(root, borderwidth=2, relief=GROOVE)
Frame4.grid(column=0, row=2, padx=25, pady=25)

## Frame 5
Frame5 = Frame(root, borderwidth=2, relief=GROOVE)
Frame5.grid(column=1, row=2, padx=25, pady=25)

## Frame 6
Frame6 = Frame(root, borderwidth=2, relief=GROOVE)
Frame6.grid(column=2, row=2, padx=25, pady=25)

## Frame 7
Frame7 = Frame(root, borderwidth=2, relief=GROOVE)
Frame7.grid(column=0, row=3, padx=25, pady=25)

## Frame 8
Frame8 = Frame(root, borderwidth=2, relief=GROOVE)
Frame8.grid(column=1, row=3, padx=25, pady=25)

## Frame 9
Frame9 = Frame(root, borderwidth=2, relief=GROOVE)
Frame9.grid(column=2, row=3, padx=25, pady=25)


Label(Frame1, text=" GENE INFORMATIONS ").pack(padx=10, pady=10)
Label(Frame2, text="BIBLIO INFORMATIONS").pack(padx=10, pady=10)
Label(Frame3, text="PROTEIN INFORMATIONS").pack(padx=10, pady=10)
Label(Frame4, text="WORK IN PROGRESS").pack(padx=10, pady=10)
#Label(Frame5, text="WORK IN PROGRESS").pack(padx=10, pady=10)
Label(Frame6, text="TEXT TO SPEAK").grid(column=1, row=0)
Label(Frame7, text="INVOLVED PATHWAYS").pack(padx=10, pady=10)
Label(Frame8, text=" GENERAL OPTIONS ").pack(padx=10, pady=10)
Label(Frame9, text="WORK IN PROGRESS").pack(padx=10, pady=10)

## Gene information Frame
## A few Labels describing the gene
label_gene_name = Label(Frame1, text="Name : NA")
label_gene_symbol = Label(Frame1, text="Symbol : NA")
label_gene_taxid = Label(Frame1, text="Taxid : NA")
label_gene_entrez = Label(Frame1, text="Entrez-Gene : NA")
label_gene_id = Label(Frame1, text="id : NA")

label_gene_name.pack()
label_gene_symbol.pack()
label_gene_taxid.pack()
label_gene_entrez.pack()
label_gene_id.pack()


## Protein Informations
## A few label describing the protein
label_uniprotId = Label(Frame3, text="UniProt ID : NA")
label_uniprotId.pack()
## Open the Uniprot entry
button_uniprot = Button(Frame3, text='Uniprot Entry', state='disabled', command=lambda x=1:open_uniprot_page())
button_uniprot.pack()
## Open the mint entry
button_mint = Button(Frame3, text='Mint Entry', state='disabled', command=lambda x=1:open_mint_page())
button_mint.pack()



## Settings for the online
## what are we looking for ? currently:
##	- Abstract -> pubmed crawling
##	- Other -> work in progress
label_abstract_found = Label(Frame2, text="Abstract Found : NA")
label_abstract_found.pack()

"""
fetch_settings = IntVar()
fetch_settings.set(0)
Radiobutton(Frame2, text="Abstract", variable=fetch_settings, value=1).pack(anchor=W)
Radiobutton(Frame2, text="Gene", variable=fetch_settings, value=2).pack(anchor=W)
Radiobutton(Frame2, text="Other", variable=fetch_settings, value=3).pack(anchor=W)

value = StringVar() 
value.set("Your Search")
entree = Entry(Frame2, textvariable=value, width=10)
button_search = Button(Frame2, text='Search', command=lambda x=2:search(value.get(), fetch_settings.get()))	
button_report = Button(Frame2, text='Report', state=DISABLED, command=lambda x=0:generate_report.write_abstract_report())

entree.pack()
button_search.pack()
button_report.pack()
"""


## Involved Pathways
## Scrollable listbox of pathways
scrollbar_y = Scrollbar(Frame7)
scrollbar_y.pack(side='right', fill='y')
scrollbar_x = Scrollbar(Frame7,  orient=HORIZONTAL)
scrollbar_x.pack(side='bottom', fill='x')
listbox_pathways = Listbox(Frame7, bg='white', yscrollcommand=scrollbar_y.set)
scrollbar_y.config(command=listbox_pathways.yview)
scrollbar_x.config(command=listbox_pathways.xview)
listbox_pathways.pack()

## Show Button
button_show = Button(Frame7, text='Show', command=lambda x=1:display_pathway(listbox_pathways.get(listbox_pathways.curselection())))	
button_show.pack()


## General Settings & Quit Button
button_settings = Button(Frame8, text='Settings', command=lambda x=1:settings.open_settings(root))
button_quit = Button(Frame8, text='Quitter', command=root.quit)
button_settings.pack()
button_quit.pack()

## Center Frame 
## include image
## Main Search -- IN PROGRESS
img = ImageTk.PhotoImage(Image.open("bibot_logo_small.png"))
panel = Label(Frame5, image = img)
panel.grid(column=1, row=0)

general_search = StringVar() 
general_search.set("WORK IN PROGRESS")
entree = Entry(Frame5, textvariable=general_search, width=18)
entree.grid(column=1, row=1)


fetch_settings_general = IntVar()
fetch_settings_general.set(0)
radiobutton_option_1 = Radiobutton(Frame5, text="Abstract", variable=fetch_settings_general, value=1)
radiobutton_option_2 = Radiobutton(Frame5, text="Gene", variable=fetch_settings_general, value=2)
radiobutton_option_3 = Radiobutton(Frame5, text="Other", variable=fetch_settings_general, value=3)

radiobutton_option_1.grid(column=0, row=2)
radiobutton_option_2.grid(column=1, row=2)
radiobutton_option_3.grid(column=2, row=2)

button_search_general = Button(Frame5, text='Search', command=lambda x=2:search(general_search.get(), fetch_settings_general.get()))	
button_report_general = Button(Frame5, text='Report', state=DISABLED, command=lambda x=0:generate_report.write_abstract_report_dev())
button_random_general = Button(Frame5, text='Reset', command=reset_search)

button_search_general.grid(column=0, row=3)
button_report_general.grid(column=1, row=3)
button_random_general.grid(column=2, row=3)



## VOICE FRAME
## IN PROGRESS
button_speak = Button(Frame6, text='Speak', command=voice.say_something)
voice_language = IntVar()
voice_language.set(1)
radiobutton_language_1 = Radiobutton(Frame6, text="English", variable=voice_language, value=1)
radiobutton_language_2 = Radiobutton(Frame6, text="French", variable=voice_language, value=2)

radiobutton_language_1.grid(column=0, row=1)
radiobutton_language_2.grid(column=2, row=1)
button_speak.grid(column=1, row=2)



root.mainloop()


