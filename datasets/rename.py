import numpy as np 
import os 

data_dir = "./boussinesq2d_vtk_v/"
files = [f for f in os.listdir(data_dir) if f.endswith(".vtk")]
for f, file in enumerate(files):
    substring = file[12: len(file) - 4]
    temp = substring.split(".")
    # print(file, temp[0], temp[1])
    newname = data_dir + "V_" + str(temp[0]).zfill(2) + "." + str(temp[1]).zfill(2) + ".vtk"
    os.rename(data_dir + file, newname)