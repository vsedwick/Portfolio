#Contact Victoria for assistance
#email: victoria.sedwick@einsteinmed.edu
#number: 317-331-1322

####Description##############
#This code will run through every folder that contains a 'SNB_Counts.csv' 
# file and compile the necessary information into an excel sheet.

#store the address for the parent folder
rootdir= r"E:\RNAScope\CRFR2-VGAT-VGLUT\Males\crfr2-vglut-vgat\Stacks\CROPPED\MeP" #UPDATE
project='RPE-test' #UPDATE

channel1='vgat'  #UPDATE
channel1_cutoff=0 #UPDATE

channel2='vglut' #UPDATE
channel2_cutoff=0 #UPDATE

channel3='crfr2' #UPDATE
channel3_cutoff=0 #UPDATE

#############Instructions###########
#This code was designed as though all stacks that go through matlab exists in an individual folder labeled by the name of the stack.
#Each folder should contain a single "..._SNB_counts.csv". All other contents will be ignored.
#This code ignores all headers. Channels must be in columns F:H in excel and Dapi must be column A
#The only information that should be changed is rootdir ln[81] and special_total ln[41-48] (if you're counting by specific criteria such as any cell >1),
#And headers in summary data frame ln[145-161]
#This code is designed to work for no more than 3 channels besides DAPI


import pandas as pd
import numpy as np
import os
from datetime import datetime


#General counts
def count_total(c,a):
    count=0
    for i in c:
        if i>a:
            count+=1
        else:
            continue
    return count

#conditional counts
#for colocalizations x is the total number in a channel, y is the channel you want to count
def count_coloc(x,y,a):
    count=0
    for j in y[0:int(x)]:
        if j>a:
            count+=1
        else:
            continue
    return count

def triple_coloc(x,y,z,b,c):
    count=0
    for i,j in zip(y[0:int(x)], z[0:int(x)]):
        if i>b and j>c:
            count+=1
        else:
            continue
    return count

#creates time stamp
def get_timestamp():
    dt= str(datetime.now())
    dt=dt.split(' ')
    dt1=dt[0]
    dt2=dt[1].split('.')
    dt2=str(dt2[0].replace(':','-'))
    ts=dt1+'_'+dt2
    return ts

#Adds an additional row for the sum of all values
def sums(a,b,c,d,e,f,g,h,i):
    a.append('SUMS')
    b.append(sum(b))
    c.append(sum(c))
    d.append(sum(d))
    e.append(sum(e))
    f.append(sum(f))
    g.append(sum(g))
    h.append(sum(h))
    i.append(sum(i))
def perc(x,y):
    try:
        z=x/y[-1]*100
    except ZeroDivisionError:
        z=0
    return z
def main():
    global rootdir
    #empty lists to be filled
    total_nuclei=[]
    names=[]
    c1_c2=[]
    c1_c3=[]
    c2_c3=[]
    c1_c2_c3=[]
    p121=[]
    p122=[]
    p131=[]
    p133=[]
    p232=[]
    p233=[]
    c1=[]
    c2=[]
    c3=[]
    
    #lists everthing inside the parent folder including individual files
    stacks=os.listdir(rootdir)
    #Interates over every folder in parent path
    for folder in stacks:
        #skips any non-folders/files with the following extensions
        if folder.endswith('.csv') or folder.endswith('.txt') or folder.endswith('.jpg') or folder.endswith('.pzfx') or folder.endswith('.xlsx') or folder.endswith('.gif'):
            continue
        else:
            #Appends folder name to path
            f=os.path.join(rootdir, folder)
            #Iterates over every file in folder
            for file in os.listdir(f):
                #file we need
                if '_SNB_Counts' in file:
                    print(file)
                    #add folder names to list
                    names.append(folder)
                    #appends file name to path
                    s=os.path.join(f, file)
                    #skip rows removes headers
                    counts=pd.read_csv(s, header=None, skiprows=[0])
                    #Total nuclei
                    total_nuclei.append(counts[0].max())
                    #count each channel
                    c1_total=(count_total(counts[5], channel1_cutoff))
                    c2_total=(count_total(counts[6], channel2_cutoff))
                    if len(counts.columns)>=8:
                        c3_total=(count_total(counts[7], channel3_cutoff))
                    #adds to list
                    c1.append(c1_total)
                    c2.append(c2_total)
                    if len(counts.columns)<8:
                        c3.append(np.nan)
                    else:
                        c3.append(c3_total)
                    #Sort sheet by c1 highest to lowest
                    up_c1=counts.sort_values(by=5, ascending=[False])
                    ##Colocalizations with C1 saved to list
                    c1c2=count_coloc(c1_total, up_c1[6], channel2_cutoff)
                    c1_c2.append(c1c2)
                    p121.append(perc(c1c2, c1))
                    p122.append(perc(c1c2, c2))
                    if len(counts.columns)<8:
                        c1_c2_c3.append(np.nan)
                    else:
                        c1_c2_c3.append(triple_coloc(c1_total, up_c1[6], up_c1[7], channel2_cutoff, channel3_cutoff))
                        c1c3=count_coloc(c1_total, up_c1[7], channel3_cutoff)
                        #colocalizations by C2 saved to list
                        up_c2=counts.sort_values(by=6, ascending=[False])
                        c2c3=count_coloc(c2_total, up_c2[7], channel3_cutoff)

                    if len(counts.columns)<8:
                        p131.append(np.nan)
                        p133.append(np.nan)
                        p232.append(np.nan)
                        p233.append(np.nan)
                    else:
                        p131.append(perc(c1c3,c1))
                        p133.append(perc(c1c3,c3))
                        p232.append(perc(c2c3,c2))
                        p233.append(perc(c2c3,c3))
                        c1_c3.append(c1c3)
                        c2_c3.append(c2c3)
                    # c1_c2_c3.append(c1c2c3)

                    #colocalizations by C2 saved to list
                    
                    #single counts
                
                else:
                    continue
    #create time stamp
    time=get_timestamp()

    #adds row that sums all values
    #names.append('SUMS')
    # sums(names, total_nuclei, c1, c2, c3, c1_c2, c1_c3, c2_c3, c1_c2_c3)
      
    # for 3 channels
    # summary1=pd.DataFrame({'name': names, 
    #                       'Total DAPI':total_nuclei, 
    #                       f'{channel1} Total':c1
    #                     #   f'{channel2} Total':c2, 
                        #   f'{channel1}+{channel2}':c1_c2, 
                        #   f'%{channel1}+{channel2}/{channel1}':p121,
                        # #   f'%{channel1}+{channel2}/{channel2}':p122, 
                        #   })

    summary1=pd.DataFrame({'name': names, 
                          'Total DAPI':total_nuclei, 
                          f'{channel1} Total':c1, 
                          f'{channel3} Total':c3, 
                          f'{channel2} Total':c2, 
                          f'{channel1}+{channel3}':c1_c3, 
                          f'{channel1}+{channel2}':c1_c2, 
                          f'{channel2}+{channel3}':c2_c3,
                          f'{channel1}+{channel2}+{channel3}':c1_c2_c3,
                          f'%{channel1}+{channel3}/{channel1}':p131,
                          f'%{channel1}+{channel3}/{channel3}':p133, 
                          f'%{channel1}+{channel2}/{channel1}':p121, 
                          f'%{channel1}+{channel2}/{channel2}':p122,
                          f'%{channel2}+{channel3}/{channel3}':p233,
                          f'%{channel2}+{channel3}/{channel2}':p232
                          })

    summary=os.path.join(rootdir, f'{project}-summary_{time}.xlsx')
    summary1.to_excel(summary, index=False, header=True)
    #create and store csv content in list
   
   

   

    
if __name__=="__main__":
    main()



    #run function




