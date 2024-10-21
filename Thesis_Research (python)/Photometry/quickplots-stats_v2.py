
summary_folder = r"C:\Users\sedwi\Desktop\New_code_tester\V. Males Non-Attacking\Summary"

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from statannotations.Annotator import Annotator
from batch_photocode_v2 import make_folder
from scipy import stats



def determine_ylabel(file_name, sheet_name):
    if "AUC Magnitude" in file_name:
        return "Δ AUC Magnitude"
    elif "AUC" in file_name:
        return "AUC"
    elif "Percent Change" in file_name:
        return "% Change"
    elif "Raw Means Differences" in file_name:
        if sheet_name in "Normalized":
            return "\Delta$F/F (event-baseline)"
        elif sheet_name in "Scaled":
            return "scaled \Delta$F/F (event-baseline)"
        elif sheet_name in "Zscore":
            return "zF (event-baseline)"
    elif "Raw Means" in file_name:
        if sheet_name in "Normalized":
            return "\Delta$F/F"
        elif sheet_name in "Scaled":
            return "scaled \Delta$F/F"
        elif sheet_name in "Zscore":
            return "zF"
    elif "Zscore" in file_name:
        return '\Delta$ Z-Score'
        
def restr_with_baseline(sheet):
    headers = list(sheet.columns)
    label = []
    behavior = []
    values = []

    for i in range(len(headers) - 1):
        if '_baseline' in headers[i]:
            print(headers[i])
            base_bar = [x for x in sheet[headers[i]] if str(x) != 'nan']
            values.extend(base_bar)
            label.extend(["Baseline" for x in base_bar])
    

            event_bar = [x for x in sheet[headers[i+1]] if str(x) != 'nan']
            values.extend(event_bar)
            label.extend(["Event" for x in event_bar])

            behavior.extend([headers[i + 1] for x in event_bar])
            behavior.extend([headers[i + 1] for x in base_bar])
    
    
    return {"Values": values, "Placement": label, "Behaviors": behavior}

def restructure_single(sheet):
    headers = list(sheet.columns)
    values = []
    behavior = []
    for i in range(len(headers) - 1):
        if 'ID' in headers[i]:
            list1 = [x for x in sheet[headers[i + 1]] if str(x) != 'nan']
            values.extend(list1)
            behavior.extend([headers[i + 1] for x in list1])

    return {"Values": values, "Behaviors": behavior}

def paired_ttest_bar(new_dict, ylabel, file_path, save_path, sht): #https://seaborn.pydata.org/tutorial/axis_grids.html 
    file_name = file_path.split("\\")[-1]
    stage = file_path.split("\\")[-3]

    df = pd.DataFrame(new_dict)
    df.head()

    args = dict(x="Placement", y="Values", order=['Baseline', 'Event'], col="Behaviors")

    # Create the grid of plots
    g = sns.catplot(data=df, sharex=False, sharey='col', edgecolor="black", 
                    errcolor="black", errwidth=1.5, capsize=0.1, 
                    palette=['darkgrey', 'deepskyblue'], height=2, 
                    col_wrap=5, aspect=1, alpha=0.3, kind="bar", 
                    ci=68, **args, margin_titles=True)

    g.map(sns.stripplot, "Placement", "Values", color='k', alpha=0.7)
    g.set_ylabels(label=ylabel)
    g.set_titles("{col_name}")
    g.figure.subplots_adjust(top=0.9)

    # Annotate each subplot with paired t-test results
    for name, ax in g.axes_dict.items():
        annot = Annotator(ax, [("Baseline", "Event")], **args, 
                          data=df.loc[df['Behaviors'] == name, :], 
                          hide_non_significant=True)
        annot.configure(test='t-test_paired').apply_and_annotate()

    # Set the figure title and save the whole grid
    name = file_name.split('.')[0]
    g.figure.suptitle(f'{name}_{sht}')

    # Save the full figure, which is the entire grid
    plt.savefig(f"{save_path}\\{name}_{sht}.tif", bbox_inches='tight')


def one_sample_Tplot(data_new_dict, ylabel, file_path, save_path, sht):

    file_name = file_path.split("\\")[-1]


    df = pd.DataFrame(data_new_dict)
    
    # plt.figure(figsize = (6, 10))
    # Create the plot
    g = sns.catplot(data=df, x="Behaviors", y="Values", kind="box", ci=68, height = 4, aspect = 2.5, showfliers = False)
    sns.stripplot(data=df, x="Behaviors", y="Values", color="k", alpha=0.5, jitter = True)
    
    # Get unique behaviors
    behav = df["Behaviors"].unique()
    
    # Get the axis object from the FacetGrid
    ax = g.ax
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Perform one-sample t-tests and annotate manually
    for i, b in enumerate(behav):
        values = df[df["Behaviors"] == b]["Values"]
        t_stat, p_value = stats.ttest_1samp(values, popmean=0)
        
        # Annotate the plot
        x_coord = i
        y_coord = values.max() + (values.max() * 0.1) 
        if p_value > 0.1 or p_value == np.NaN:
            annotation = ''
        else:
            annotation = f"p = {p_value:.3f}"
        ax.text(x_coord, y_coord, annotation, ha='center')
    
    # Set labels and titles
    ax.set_ylabel(ylabel)
    name = file_name.split('.')[0]
    ax.set_title(f"{name}_{sht}")
    
    # Save the figure
    g.savefig(f"{save_path}/{name}_{sht}.tif")
    
def main(summary_folder):
    stages = os.listdir(summary_folder)
    for stage in stages:
        stage_dir = os.path.join(summary_folder, stage)
        resolution = os.listdir(stage_dir)
        for i in resolution:
            if i not in 'Figures' or i not in 'zscore':
                test_dir = os.path.join(stage_dir, i)
                to_analyze = os.listdir(test_dir)
                for file_name in to_analyze:
                    if file_name not in "Figures":
                        file_path = os.path.join(test_dir, file_name)
                        wb = pd.ExcelFile(file_path)
                        for sht in wb.sheet_names:
                            sheet = pd.read_excel(file_path, sheet_name= sht)
                            # print(sheet)
                            ylabel = determine_ylabel(file_name,sht)
                            figure_folder = make_folder("Figures", test_dir)
                            save_path = make_folder("Prelim bar plots", figure_folder)
                            check_if_baseline = [i for i in sheet.columns if '_baseline' in i]
                            if len(check_if_baseline) == 0:
                                new_dict = restructure_single(sheet)
                                one_sample_Tplot(new_dict, ylabel, file_path, save_path, sht)
                            else:
                                print("Plotting double bar")
                                new_dict = restr_with_baseline(sheet)
                                print(len(new_dict["Values"]), len(new_dict["Placement"]), len(new_dict["Behaviors"]))
                                if stage == 'during':
                                    df = pd.DataFrame(new_dict)
                                    df.to_csv(f'{figure_folder}\\new_dict.csv')
                                paired_ttest_bar(new_dict, ylabel, file_path, save_path, sht)

if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(summary_folder)