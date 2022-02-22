#Fixes paths in DFN input files
#adapted from DFNWorks fix_paths.py

import os

old = 'FILEPATH'
new = os.path.join('/',*os.getcwd().split('/'))
path = os.getcwd()
files = os.listdir(path)

for file in files:
    if file == "fix_paths.py":
        continue
    try:
        with open(file,'r') as f:
            replaced_text = f.read().replace(old,new)
    except UnicodeDecodeError:
        continue
    with open(file,'w') as f:
        f.write(replaced_text)
