"""
Setting stuff
"""
from Tkinter import *

def open_settings(root):
	"""
	-> Open the settings window
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
	button_quit = Button(SettingsFrame2, text='Quitter', command=settings_window.destroy)
	button_save.pack()
	button_default.pack()
	button_quit.pack()
	SettingsFrame2.pack(side=LEFT, padx=30, pady=30)


