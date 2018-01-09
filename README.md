# BIBOT
a few things for biblio

http://www.nltk.org/book/ch07.html


# Main Project
=> Parse input
	- detect if it's more loikely to be a gene, uniprot id, etc.
=> Screen a few database
	- pubmed (fetch abstract & articles)
	- mint (molecular interaction db)
	- kegg (metabolic pathway)
=> Create a report for the user
	- html
	- pdf (latex)
=> GUI
	- use a tkinter graphical user interface

# Additionnal stuff

=> Reasonner
-Problème du sens des phrases, comprendre les relations entre les sujets.
	-> Solution 1: liste exhaustives des relations

-Tester pertinence d'un ensemble de relation

=> Parser (ToDo)
-Taguer les sujets
-Taguer les relations

=> Graphe representation
	- deal with sif & GDF file
	- draw interactions

=> Text to speach 

# Plan

1) Retrieve Article (Absrtract) from NCBI
 -> Done with the fetch_abstract function in trashlib
2) Parse Abstract
 -> Working on it with nltk
3) Reason on parsed abstract

## Requirements:
Biopython
nltk
mygene 3.0.0 (python module)
PyEntrezId 1.5.8
Pillow
pyttsx (32 bit version for windows)
Microsoft Visual C++ 9 under windows

## Installation:
pip install biopython
pip install bioservices
pip install nltk
pip install mygene
pip install pyttsx
