#!/usr/bin/python3

import json

f = open('configs/vertexSelection_2016_simp_preselection.json')

data = json.load(f)

for i, (cut,info) in enumerate(data.items()):
    print(cut)
    print(info)
    newdict = {}
    cut_id  = info['id']
    for othcut, othinfo in data.items():
        if othinfo['id'] <= cut_id:
            newdict[othcut] = othinfo
    with open('selections/cutflow/%s_selection_inclusive.json'%(cut), "w") as outfile:
        json.dump(newdict, outfile, indent = 4)
    

