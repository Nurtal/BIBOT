"""
Graphical Interface for BIBOT
"""


from Tkinter import *
from bibliosearch import *



"""Test Stuff"""

def testFunction(machin):
	print "Testing "+str(machin)




"""Interface"""

root = Tk()
root.title('BIBOT Interface')
root['bg']='white'





# frame 1
Frame1 = Frame(root, borderwidth=2, relief=GROOVE)
Frame1.pack(side=LEFT, padx=30, pady=30)

# frame 2
Frame2 = Frame(root, borderwidth=2, relief=GROOVE)
Frame2.pack(side=LEFT, padx=10, pady=10)

# frame 2
Frame3 = Frame(root, borderwidth=2, relief=GROOVE)
Frame3.pack(side=RIGHT, padx=30, pady=30)


Label(Frame1, text="Frame 1").pack(padx=10, pady=10)
Label(Frame2, text="Frame 2").pack(padx=10, pady=10)
Label(Frame3, text="Frame 3").pack(padx=10, pady=10)

value = StringVar() 
value.set("Your Search")
entree = Entry(Frame2, textvariable=value, width=10)
searchButton = Button(Frame2, text='Search', command=lambda x=1:fetch_abstract(value.get()))	





qb = Button(Frame3, text='Quitter', command=root.quit)


entree.pack()
searchButton.pack()
qb.pack()

root.mainloop()


