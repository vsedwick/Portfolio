import pandas as pd
import os
import shutil

def main():
    #navigate to parent folder
    rootdir= "F:\PGR-CRFR2\stacks"
    #lists everthing inside the directory
    stacks=os.listdir(rootdir)
    
    #empty lists to be filled
    

    #Interates over every folder in root path
    for file in stacks:
        #skips any non-folders/files with the following extensions
        if not file.endswith('.tif'):
            continue
        else:
            file_name=[]
            #source path
            source=os.path.join(rootdir, file)
            #Removes extention from file name
            f=file.replace('.tif','')
            #remove other stuff
            # f=f.split('-')
            file_name.append(f)
            # file_name.append(f[1])
            folder_name='-'.join(file_name)

            #Creates path with file name
            folder_source=os.path.join(rootdir, folder_name)
            #mode (from internet https://www.geeksforgeeks.org/create-a-directory-in-python/)
            mode = 0o666
            os.mkdir(folder_source, mode)


            #move file to folder
            new_home=os.path.join(folder_source, file)
            shutil.move(source, new_home)
            

    
if __name__=="__main__":
    main()
