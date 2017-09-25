"""
Setting stuff
"""
from Tkinter import *

def open_settings(root):
	"""
	-> Open the settings window
	
	[IN PROGRESS]
	"""

	## General
	settings_window = Toplevel(root)
	settings_window.title('BIBOT Settings')

	## Frame 1 - NCBI search settings
	SettingsFrame1 = Frame(settings_window, borderwidth=2, relief=GROOVE)
	Label(SettingsFrame1, text="NCBI").pack(padx=10, pady=10)
	max_abstract_return = IntVar()
	max_abstract_return.set(5)
	Label(SettingsFrame1, text="Maximum number of \n abstracts to return").pack()
	entree = Entry(SettingsFrame1, textvariable=max_abstract_return, width=10)
	entree.pack()
	SettingsFrame1.pack(side=LEFT, padx=30, pady=30)

	## Frame 2 - Save & Quit
	SettingsFrame2 = Frame(settings_window, borderwidth=2, relief=GROOVE)
	label_status = Label(SettingsFrame2, text="Changes Unsaved").pack(padx=10, pady=10)
	button_save = Button(SettingsFrame2, text='Save', command=lambda x=2:write_settings("custom", max_abstract=max_abstract_return.get()))
	button_default = Button(SettingsFrame2, text='Default', command=lambda x=1:write_settings("default"))
	button_quit = Button(SettingsFrame2, text='Quitter', command=settings_window.destroy)
	button_save.pack()
	button_default.pack()
	button_quit.pack()
	SettingsFrame2.pack(side=LEFT, padx=30, pady=30)



def write_settings(mode, **custom_values):
	"""
	-> Write the settings in a settings.csv file
	-> mode is a string:
		- default : use default values
		- custom : use values contain in **custom_values
	-> custom values : set of keyword arguments:
		-max_abstract : the max number of abstract to return with a classic pubmed search
	
	[IN PROGRESS]
	"""

	## Init default values
	max_abstract_return = 5

	## Custom case
	if(mode == "custom"):

		## Get optionnal parameter (custom values)
		if('max_abstract' in custom_values):
			max_abstract_return = custom_values['max_abstract']

	## Write settings file
	settings_file = open("settings.csv", "w")
	settings_file.write("max_abstract_return,"+str(max_abstract_return)+"\n")
	settings_file.close()

	## Update the label status
	## [TODO]


def read_settings():
	"""
	-> Read the settings file (settings.csv) and return
	   a dictionnary where key are the parmater names
	   and values the values of parameters
	"""

	## init dictionnary to return
	settings_to_values = {}

	## read the setting file
	settings_file = open("settings.csv", "r")
	for line in settings_file:
		line = line.replace("\n", "")
		line_in_array = line.split(",")
		settings_to_values[line_in_array[0]] = line_in_array[1]
	settings_file.close()

	## return dictionnary
	return settings_to_values

