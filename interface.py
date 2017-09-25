"""
Graphical Interface for BIBOT
"""


from Tkinter import *
#import ImageTk
from bibliosearch import *
import input_parser
import generate_report
import settings


"""Test Stuff"""

def testFunction(machin):
	print "Testing "+str(machin)


def search(query, query_type):
	"""
	Wrapper around the search parser
	from input_parser
	"""
	input_parser.search_parser(query, query_type)
	button_report['state'] = "normal"

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
			if(line_in_array[0] == "uniprot_id"):
				label_uniprotId["text"] = "UniProt ID : "+str(line_in_array[1])
		protein_data_file.close()


"""Interface"""

root = Tk()
#img = ImageTk.PhotoImage(Image.open("bibot_logo.png"))
root.title('BIBOT Interface')
root['bg']='white'


## main frame 

# frame 1
Frame1 = Frame(root, borderwidth=2, relief=GROOVE)
Frame1.grid(column=0, row=1)

# frame 2
Frame2 = Frame(root, borderwidth=2, relief=GROOVE)
Frame2.grid(column=1, row=1)

# frame 3
Frame3 = Frame(root, borderwidth=2, relief=GROOVE)
Frame3.grid(column=2, row=1)

## Frame 4
Frame4 = Frame(root, borderwidth=2, relief=GROOVE)
Frame4.grid(column=0, row=2)

## Frame 5
Frame5 = Frame(root, borderwidth=2, relief=GROOVE)
Frame5.grid(column=1, row=2)

## Frame 6
Frame6 = Frame(root, borderwidth=2, relief=GROOVE)
Frame6.grid(column=2, row=2)

Label(Frame1, text="Gene Information").pack(padx=10, pady=10)
Label(Frame2, text="Biblio Search").pack(padx=10, pady=10)
Label(Frame3, text="Protein Information").pack(padx=10, pady=10)
Label(Frame4, text="Involved Pathways").pack(padx=10, pady=10)
Label(Frame5, text="General Options").pack(padx=10, pady=10)
Label(Frame6, text="WORK IN PROGRESS").pack(padx=10, pady=10)


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


## Settings for the online
## what are we looking for ? currently:
##	- Abstract -> pubmed crawling
##	- Other -> work in progress
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

## Involved Pathways
## IN PROGRESS
scrollbar_y = Scrollbar(Frame4)
scrollbar_y.pack(side='right', fill='y')

scrollbar_x = Scrollbar(Frame4,  orient=HORIZONTAL)
scrollbar_x.pack(side='bottom', fill='x')


listbox_pathways = Listbox(Frame4, bg='white', yscrollcommand=scrollbar_y.set)
scrollbar_y.config(command=listbox_pathways.yview)
scrollbar_x.config(command=listbox_pathways.xview)

#scrollbar = Scrollbar(listbox_pathways, orient=VERTICAL)
#scrollbar.pack()
listbox_pathways.pack()










## General Settings & Quit Button
button_settings = Button(Frame5, text='Settings', command=lambda x=1:settings.open_settings(root))
button_quit = Button(Frame5, text='Quitter', command=root.quit)

button_settings.pack()
button_quit.pack()

root.mainloop()


