import sys
import pandas as pd
import os
import numpy as np
import openpyxl
from combined_photocode_02262024 import open_gui


behavior_row = 36
lolo = 'start'
fps=30
length_s=300

def main():

    global behavior_row, lolo, fps, length_s;

    rootdir= r"D:\Experiments\Lesion\Cohort 1-(date lesion)\Pup Exposure\Lesion Pup Exposure\Export Files"

    #lists everthing inside the parent folder including individual files
    ext=[".avi",".csv",".mp4", ".pzfx", ".jpg", ".txt", ".prism", ".xlsx", ".doc"]
    excluded_folders = ["Behaviors_PA", "Behaviors", "Summary", "Summary_PA", "Archive", "Videos", 'Behavior_StatsSummary.xlsx']
    folder_list = [
        folder for folder in os.listdir(rootdir)
        if folder not in excluded_folders or folder.lower().endswith(tuple(ext))
    ]
    print(folder_list)
    #Interates over every folder in parent path
    behavior = []
    cum = {}
    freq = {}
    lat = {}
    lat_last = {}
    m = 1

    
    for folder in folder_list:
        # print(folder)
        folder_path = os.path.join(rootdir, folder)
        files = os.listdir(folder_path)
        for file in files:
            if '$' in file:
                continue
            elif 'Raw data' in file and file.endswith('.xlsx') and '$' not in file:
                s=os.path.join(folder_path,file)

                print(s)
            ##Read sheet for behavior
                animal=pd.read_excel(s, sheet_name = [0], header = [behavior_row-1], skiprows=[behavior_row])
                animal = animal[0]

                behaviors = animal.columns
                if lolo not in behaviors:
                    # new_row=int(input("Behaviors are not in designated row. Try again: "))
                    animal = pd.read_excel(s, header=[33-1], skiprows=[33])

                while m<=1:
                    behaviors = open_gui(file, behavior_row)
                    behavior.extend(behaviors)
                    m+=1
                start,end = crop(animal, fps, length_s, lolo)
                animal = animal[start:end]
                
                # ##Calculations.
                cum1 = [cumulative(animal[b])/fps for b in behavior]
                lat1 = [latency(np.array(animal[b]))/fps for b in behavior]
                freq1 = [frequency(np.array(animal[b])) for b in behavior]
                lat_last1 = [latency_last(np.array(animal[b]))/fps for b in behavior]
                
                cum.update({f"{folder}": cum1})
                lat.update({f"{folder}": lat1})
                freq.update({f"{folder}": freq1})
                lat_last.update({f"{folder}": lat_last1})
    

    cummulative1 = pd.DataFrame(cum).T
    cummulative1.columns = behavior

    latency1 = pd.DataFrame(lat).T
    latency1.columns = behavior

    frequency1 = pd.DataFrame(freq).T
    frequency1.columns = behavior

    latency_last1 = pd.DataFrame(lat_last).T
    latency_last1.columns = behavior

    with pd.ExcelWriter(f'{rootdir}\\StatsSummary.xlsx', engine = "openpyxl") as writer:
        cummulative1.to_excel(writer, sheet_name = "Cumulative Duration")

        latency1.to_excel(writer, sheet_name = "Latency")

        frequency1.to_excel(writer, sheet_name = "Frequency")

        latency_last1.to_excel(writer, sheet_name = "Latency to Last")
        # writer.save()
        # writer.close()

def crop(animal, fps, length_s, lolo):
    
    start = 0
    for i in range(len(animal[lolo])):
        if animal[lolo][i] == 1:
            start+=i
            break
    end = start + (length_s*fps) + 1
    return start, end;

def latency(x):
    m = 0
    if 1 not in x:
        return length_s*fps
    for i in range(len(x)):
        if x[i] == 1:
            m+=i
            return m
        
def latency_last(x):
    ones_place = [i for i in range(len(x)) if x[i] == 1]
    if len(ones_place) > 0:
        return ones_place[-1]
    else:
        return 0
        
def cumulative(x):
    k=0
    for i in x:
        if i==1:
            k+=1

    return k

def frequency(x):
    k=0
    a=x[:-1]
    b=x[1:]
    for i,j in zip(a,b):
        if j==1 and i==0:
            k+=1
        else:
            continue
    return k

def make_folder(x, project_home):

    mode = 0o666
    j=os.listdir(project_home)
    if x not in j:
        behavior_path=os.path.join(project_home, x)
        os.mkdir(behavior_path, mode)
    else:
        behavior_paths=os.path.join(project_home, x)
    behavior_paths=os.path.join(project_home, x)
    return behavior_paths



if __name__=="__main__":
    main()