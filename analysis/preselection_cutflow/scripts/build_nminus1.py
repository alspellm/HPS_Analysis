#!/usr/bin/python3

import json

f = open('configs/vertexSelection_2016_simp_preselection.json')

data = json.load(f)

for cut,info in data.items():
    print(cut)
    print(info)
    newdict = {}
    cut_id  = 0
    for othcut, othinfo in data.items():
        if othcut == cut:
            continue
        othinfo['id'] = cut_id
        newdict[othcut] = othinfo
        cut_id = cut_id + 1
    with open('selections/%s_selection.json'%(cut), "w") as outfile:
        json.dump(newdict, outfile, indent = 4)
        print(newdict)
    

