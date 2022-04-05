import numpy as np 
import os 

data_dir = "./tangaroa3d_vtk/W/"
files = [f for f in os.listdir(data_dir) if f.endswith(".vtk")]
for f, file in enumerate(files):
    substring = file[10: len(file) - 4]
    temp = substring.split(".")
    # print(file, temp[0], temp[1])
    newname = data_dir + "W_" + str(temp[0]).zfill(2) + "." + str(temp[1]).zfill(2) + ".vtk"
    print(file, newname)
    os.rename(data_dir + file, newname)