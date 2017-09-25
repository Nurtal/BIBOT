"""
Graphical Interface for BIBOT
"""


from Tkinter import *
#import ImageTk
from bibliosearch import *
import input_parser
import generate_report


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



def open_settings():
	"""
	Work in Progress
	"""

	## General 
	settings_window = Toplevel(root)
	settings_window.title('BIBOT Settings')

	## Frame 1 - NCBI search settings
	SettingsFrame1 = Frame(settings_window, borderwidth=2, relief=GROOVE)
	Label(SettingsFrame1, text="NCBI").pack(padx=10, pady=10)
	value = IntVar()
	value.set(5)
	Label(SettingsFrame1, text="Maximum number of \n abstracts to return").pack()
	entree = Entry(SettingsFrame1, textvariable=value, width=10)
	entree.pack()
	SettingsFrame1.pack(side=LEFT, padx=30, pady=30)

	## Frame 2 - Save & Quit
	SettingsFrame2 = Frame(settings_window, borderwidth=2, relief=GROOVE)
	Label(SettingsFrame2, text="Save & Quit").pack(padx=10, pady=10)
	button_save = Button(SettingsFrame2, text='Save', command=settings_window.quit)
	button_default = Button(SettingsFrame2, text='Default', command=settings_window.quit)
	button_quit = Button(SettingsFrame2, text='Quitter', command=settings_window.quit)

	button_save.pack()
	button_default.pack()
	button_quit.pack()

	SettingsFrame2.pack(side=LEFT, padx=30, pady=30)



"""Interface"""

root = Tk()
#img = ImageTk.PhotoImage(Image.open("bibot_logo.png"))

root.title('BIBOT Interface')
root['bg']='white'



# frame 1
Frame1 = Frame(root, borderwidth=2, relief=GROOVE)
Frame1.pack(side=LEFT, padx=30, pady=30)

# frame 2
Frame2 = Frame(root, borderwidth=2, relief=GROOVE)
Frame2.pack(side=LEFT, padx=10, pady=10)

# frame 3
Frame3 = Frame(root, borderwidth=2, relief=GROOVE)
Frame3.pack(side=RIGHT, padx=30, pady=30)


Label(Frame1, text="Gene Info").pack(padx=10, pady=10)
Label(Frame2, text="Biblio Search").pack(padx=10, pady=10)
Label(Frame3, text="Options").pack(padx=10, pady=10)


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


## General Settings & Quit Button
button_settings = Button(Frame3, text='Settings', command=open_settings)
button_quit = Button(Frame3, text='Quitter', command=root.quit)

button_settings.pack()
button_quit.pack()

root.mainloop()


