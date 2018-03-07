
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


import numpy as np

## get data
category_list = ["Diagnostic", "Therapeutic", "Modelisation", "Unclassified"]
category_to_tech = {}
data = open("rush_data.txt", "r")
for line in data:
    line = line.replace("\n", "")
    line_in_array = line.split(",")
    
    if(line_in_array[0] not in category_to_tech.keys()):
        category_to_tech[line_in_array[0]] = {}

    category_to_tech[line_in_array[0]][line_in_array[1]] = int(line_in_array[2])

data.close()

## get all keys for all category
tech_list = []
for category in category_to_tech.keys():
    for tech in category_to_tech[category].keys():
        if(tech not in tech_list):
            tech_list.append(tech)

for category in category_to_tech.keys():
    for tech in tech_list:
        if(tech not in category_to_tech[category].keys()):
            category_to_tech[category][tech] = 0


print category_to_tech

## deal with low values
tech_to_remove = []
for tech in tech_list:
    d_tech_score = category_to_tech["Diagnostic"][tech]
    t_tech_score = category_to_tech["Therapeutic"][tech]
    m_tech_score = category_to_tech["Modelisation"][tech]

    if(d_tech_score < 2 and t_tech_score < 2 and m_tech_score < 2):
        tech_to_remove.append(tech)


for tech in tech_to_remove:
    del category_to_tech["Diagnostic"][tech]
    del category_to_tech["Therapeutic"][tech]
    del category_to_tech["Modelisation"][tech]

x = category_to_tech["Diagnostic"].keys()
d = category_to_tech["Diagnostic"].values()
t = category_to_tech["Therapeutic"].values()
m = category_to_tech["Modelisation"].values()

## Figure 4
N = len(x)
ind = np.arange(N)  # the x locations for the groups
ind = np.arange(0, N * 10, 10)
width = 1.5       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)

yvals = category_to_tech["Diagnostic"].values()
rects1 = ax.bar(ind, yvals, width, color='r')
zvals = category_to_tech["Therapeutic"].values()
rects2 = ax.bar(ind+width, zvals, width, color='g')
kvals = category_to_tech["Modelisation"].values()
rects3 = ax.bar(ind+width*2, kvals, width, color='b')
ax.set_ylabel('Articles')
ax.set_xticks(ind+width)
ax.set_xticklabels( category_to_tech["Diagnostic"].keys(), fontsize=13)
plt.xticks(rotation=75)
ax.legend( (rects1[0], rects2[0], rects3[0]), ('Diagnostic', 'Therapeutic', 'Modelisation') )

plt.close()



## Figure 3

category_to_tech = {}

category_to_tech["Diagnostic"] = {}
category_to_tech["Therapeutic"] = {}
category_to_tech["Modelisation"] = {}
category_to_tech["All"] = {}
category_to_tech["Diagnostic"]["regression"] = 38.9
category_to_tech["Diagnostic"]["machine"] = 33.3
category_to_tech["Diagnostic"]["other"] = 27.8
category_to_tech["Therapeutic"]["regression"] = 40.0
category_to_tech["Therapeutic"]["machine"] = 20.0
category_to_tech["Therapeutic"]["other"] = 40.0
category_to_tech["Modelisation"]["regression"] = 37.8
category_to_tech["Modelisation"]["machine"] = 42.2
category_to_tech["Modelisation"]["other"] = 20.0
category_to_tech["All"]["regression"] = 27.0
category_to_tech["All"]["machine"] = 19.2
category_to_tech["All"]["other"] = 53.8


for category in category_to_tech.keys():

    print category

    regression_techniques = category_to_tech[category]["regression"]
    machine_learning_techniques = category_to_tech[category]["machine"]
    other_techniques = category_to_tech[category]["other"]

    x_vector = [regression_techniques, machine_learning_techniques, other_techniques]
    labels = ["total regression techniques", "total machine learning techniques", "total other techniques"]
    
    plt.pie(x_vector, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    #plt.show()
    plt.close()

    patches, texts, autotexts = plt.pie(x_vector, labels=labels,
                                    autopct='%1.1f%%',
                                    shadow=True, startangle=90)
    plt.axis('equal')

    # Make the labels on the small plot easier to read.
    for t in texts:
        t.set_size('xx-large')
    for t in autotexts:
        t.set_size('smaller')
    autotexts[0].set_color('y')

    plt.show()