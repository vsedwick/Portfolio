import pandas as pd
import numpy as np
import os
import csv
from scipy import stats
import matplotlib.pyplot as plt
from statistics import mean
import openpyxl

compile_directory = r"E:\Photometry-Fall2022\Final Analysis\REANALYSIS_zscore_nolimit\Virgin Females"
# subjects = ['Virgin Males', 'Virgin Females', 'Mated Males', 'Mated Females'] #As written in folders
# subjects = ['FemaleBlood', 'MaleBlood', 'Pup Blood', 'Pups'] #As written in folders
# subjects = ['ball', 'pups', 'fake']
subjects = ['Juv fem', 'Castrated', 'Fake Pup', 'Pup Odor', 'Restricted Pup', 'Pup']

# subjects = ['Castrated', 'Fake Pup', 'Juv fem', 'Pup', 'Restricted Pup', 'Pup Odor']
Behaviors = ['Attack', 'Sniff', 'Aggressive Behaviors', 'Approach', 'Dead_Pup', 'Grooming', 'Aggressive groom', 'Digging', 'Rearing', 'No Stim Interaction', 'Nudge', 'Carrying', 'Stim Contact', 'Stim Interaction', 'qtip', 'stick', 'Cotton tip', 'Sniff_inZone', 'Grooming_inZone', 'Sniff_OutZone', 'In nest', 'Qtip Interaction']

# states = ['control', 'during', 'before']

# use_states = True

def main():
    rootdir = os.listdir(compile_directory)

    compile_list = []
    for i in rootdir:
        if i in subjects:
            print(i)
            # summary_folder = os.path.join(compile_directory, f"{i}\\Summary") 
            tag_list = compile_directory.split('\\')
            tag_list2 = tag_list[-1].split('_')
            tag = tag_list2[-1]
            print(tag)
            summary_folder = os.path.join(compile_directory, f'{i}\\Summary\\Behavior_analysis\\during\\Individual Values')
            summaries = os.listdir(summary_folder)
            for summary in summaries:
                summary_path = os.path.join(summary_folder, summary)
                compile_file_path = os.path.join(compile_directory, f'IndV_ID_{tag}_during_{summary}_compiled.xlsx')
                if not os.path.exists(compile_file_path):
                    create_excel(subjects, Behaviors, summary, compile_file_path)
                compile_list.append(compile_file_path)
                write_data_to_workbook(summary_path, compile_file_path, i)     
    # for c in compile_list: delete_empty_sheets(c)
    for compile_file_path in compile_list:
        delete_empty_sheets(compile_file_path)

def create_excel(subjects, behaviors, summary, compile_file_path):
    # Create a new Excel workbook
    wb = openpyxl.Workbook()
    
    # Remove default sheet
    default_sheet = wb.active
    wb.remove(default_sheet)
    
    # Create a worksheet for each behavior
    for behavior in behaviors:
        ws = wb.create_sheet(behavior)
        count = 0
        # Write subjects as headers
        for col, subject in enumerate(subjects, start=1):
            ws.cell(row = 1, column = col+count, value = f'ID_{subject}')
            count+=1
            ws.cell(row=1, column=col+count, value=subject)
        
    # Save the workbook
    wb.save(compile_file_path)
    # wb.close()


def write_data_to_workbook(csv_filepath, excel_filename, subject_name):  
    # Open the CSV file and read the header to get the behavior name
    df = pd.read_csv(csv_filepath)
    # Load the existing Excel workbook
    wb = openpyxl.load_workbook(excel_filename)
    df_columns = df.columns
    for column in df_columns:
        if column in wb.sheetnames:
            ws = wb[column]
            subject_col = None
            for col_index in range(1, ws.max_column + 1):
                if ws.cell(row=1, column=col_index).value == subject_name:
                    subject_col = col_index
                    break
            data = df[column]
            data_id = df[f'ANIMAL ID_{column}']
            for row_index, value in enumerate(data, start=2):
                ws.cell(row=row_index, column=subject_col, value=value)
            for row_index, value in enumerate(data_id, start=2):
                ws.cell(row=row_index, column=subject_col-1, value=value)
    print(f'{subject_name} Data is written to {csv_filepath}')
    # Save the workbook
    wb.save(excel_filename)



def delete_empty_sheets(filename):
    # Load the workbook
    wb = openpyxl.load_workbook(filename)
    
    # Iterate through all sheets in the workbook
    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]
        
        # Check if the sheet is empty
        if not any(ws.iter_rows(values_only=True)):
            # Delete the sheet
            wb.remove(ws)
    
    # Save the workbook
    wb.save(filename)

main()