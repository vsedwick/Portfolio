# AUTHOR: VICTORIA SEDWICKts
# ADAPTED FROM NEUROPHOTOMETRICS AND ILARIA CARTA (MATLAB)

# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
# https://stacbehav_filesoverflow.com/questions/11190735/python-matplotlib-superimpose-scatter-plots

project_home="F:\Photometry-Fall2022\Pup Block\Virgin\Trial 3-PE\Pup\Male"
# PROCEED WITH CAUTION CHANGING ANYTHING AFTER THIS LINE
# MAKE SURE ALL PACKAGES ARE INSTALLED
# from batch_peakanalysis_nocontrols import project_home
import pandas as pd
import numpy as np
import os
import csv
from scipy import stats
import matplotlib.pyplot as plt
from statistics import mean
import seaborn as sns
import scipy.stats
from statannotations.Annotator import Annotator
# import statannot

root = os.path.join(project_home, "Summary_PA")
states = os.listdir(root)

# CREATES SUMMARY FOLDER


def main():
    # ITERATES THROUGH 'DURING STIM', BEFORE STIM', 'CONTROL', ETC. FOLDERS
    for state in states:
        entrypoint = os.path.join(root, state)
        mode = 0o666
        j = os.listdir(entrypoint)
        if "Figures" not in j:
            summary_path = os.path.join(entrypoint, "Figures")
            os.mkdir(summary_path, mode)
        summary_path = os.path.join(entrypoint, "Figures")

        folders = os.listdir(entrypoint)
        for folder in folders:
            if 'zscore' in folder or  'Figures' in folder:
                continue
            else:
                files= os.path.join(entrypoint, folder)
                to_analyze=os.listdir(files)
                for analysis in to_analyze:
                    t=rf'{files}\{analysis}'
                    sheet=pd.read_csv(t, header=[0])

                    if 'DeltaF' not in analysis:
                        print(analysis)
                        headers=[i for i in sheet.columns.values]
                        Placement=[]
                        Values=[]
                        pval=[]
                        Behaviors=[]
                        for i in range(len(headers)-1):
                            if 'baseline' in headers[i]:
                                a=[x for x in sheet[headers[i]] if str(x)!='nan']; Values.extend(a) #baselines
                                Placement.extend(["Baseline" for x in sheet[headers[i]] if str(x)!='nan'])

                                b=[x for x in sheet[headers[i+1]] if str(x)!='nan'] ; Values.extend(b) #event
                                Placement.extend(["Event" for x in sheet[headers[i+1]] if str(x)!='nan'])
                                Behaviors.extend([headers[i+1] for x in sheet[headers[i]] if str(x)!='nan']) #for p in range(len(Baselines))])
                                Behaviors.extend([headers[i+1]for x in sheet[headers[i]] if str(x)!='nan'])
                        
                        val_dic = {"F/F": Values, "Placement": Placement, "Behaviors": Behaviors}
                        if 'AUC' in analysis:
                            barplot(val_dic, "AUC", analysis, summary_path)
                        if 'Percent' in analysis:
                            barplot(val_dic, r'%$\Delta$F/F', analysis, summary_path)
                        if 'Raw Means' in analysis:
                            barplot(val_dic, "F/F", analysis, summary_path)
                        # if 'Amplitude' in analysis:
                        #     barplot(val_dic, 'Amplitude', analysis, summary_path)
                    elif "DeltaF" in analysis:
                        print(analysis)
                        headers=[i for i in sheet.columns.values]
                        Values=[]
                        pval=[]
                        Behaviors=[]
                        for i in range(len(headers)-1):
                            if 'ANIMAL' in headers[i]:
                                b=[x for x in sheet[headers[i+1]] if str(x)!='nan'] ; Values.extend(b) #event

                                Behaviors.extend([headers[i+1] for x in sheet[headers[i+1]] if str(x)!='nan']) #for p in range(len(Baselines))])
                        val_dic = {"F/F": Values, "Behaviors": Behaviors}
                        plot2(val_dic, name, analysis, summary_path)

def barplot(val_dic, x, analysis, summary_path): #https://seaborn.pydata.org/tutorial/axis_grids.html

    df=pd.DataFrame(val_dic)
    df.head()
    print(df)

    args = dict(x="Placement", y="F/F", order=['Baseline', 'Event'])

    g=sns.catplot(data =df, sharex=False, sharey='col', edgecolor="black", errcolor="black", errwidth=1.5, capsize = 0.1, palette=['darkgrey', 'deepskyblue'], height=2, col_wrap=5, aspect=.7,alpha=0.3, kind="bar", ci = 68, col="Behaviors", **args)
    g.map(sns.stripplot, "Placement", "F/F", color='k', alpha=0.7)
    g.set_ylabels(label=x)
    g.set_titles("{col_name}")
    g.fig.subplots_adjust(top=0.9)

    for name, ax in g.axes_dict.items(): #https://www.appsloveworld.com/python/1054/how-does-one-insert-statistical-annotations-e-g-p-values-into-a-seaborn-figure
        # for ax in ax_n:
        annot = Annotator(ax, [("Baseline", "Event")], **args, data=df.loc[df['Behaviors']==name,:], hide_non_significant=True)
        annot.configure(test='t-test', test_args={'mu': 0}).apply_and_annotate()

    name, _=analysis.split('.')
    g.fig.suptitle(f'{name}')
    plt.show()
    g.savefig(f"{summary_path}\{name}.svg")

def plot2(val_dic, name, analysis, summary_path):
    df=pd.DataFrame(val_dic)
    print(df)
    # df = df.explode('F/F')
    df.head()
    g=sns.catplot(data=df, x="Behaviors", y="F/F",kind="box", ci = 68,);
    g=sns.stripplot(data=df, x="Behaviors", y="F/F")
    behav=[*set(val_dic["Behaviors"])]
    pairs = []


    for index in range(len(behav)) :
        for index_2 in range(len(behav)) :
            if index != index_2:
                pairs.append(tuple((behav[index], behav[index_2])))
    print(pairs)
    annot = Annotator(g, pairs, data=df, x="Behaviors", y="F/F",kind="box", hide_non_significant=True)
    annot.configure(test='t-test_welch').apply_and_annotate() #should be one sample t test
    # g.set_ylabels(label=r'$\Delta$F/F')
    # g.set_titles("All Behaviors")
    plt.show()

    name, _=analysis.split('.')
    g.savefig(f"{summary_path}\{name}.svg")
main()

        # CREATES SUB-FOLDER FOR EACH STATE IN SUMMARY FOLDER
