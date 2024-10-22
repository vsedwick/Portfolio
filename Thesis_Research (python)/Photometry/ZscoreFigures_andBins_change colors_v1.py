compile_directory = "E:\REANALYSIS_nolimit\Re-Analysis_Pups"

colors = {'Infanticidal Males': 'orange',
          'Virgin Males': "blue",
          'Virgin Females': 'magenta',
          "Mated Males": 'lime',
          "Mated Females": 'darkviolet'}

combos = [['Infanticidal Males', 'Virgin Males', 'Virgin Females', 'Mated Males', 'Mated Females'],
          ['Infanticidal Males', 'Virgin Males', 'Virgin Females'], 
          ['Virgin Females', 'Mated Females'], 
          ['Infanticidal Males', 'Virgin Males', 'Mated Males'], 
          ['Infanticidal Males', 'Virgin Males'], 
          ['Virgin Males', 'Mated Males'],
          ['Infanticidal Males','Mated Males'],
          ['Infanticidal Males','Mated Females'],
          ['Mated Males', 'Mated Females']]

project1 = ['Pup Exposure', 'Virgins', 'Females', 'Males', 'Virgin Males', 'Non Aggressive Males', 'VMI v MM', 'VMI v MF', 'Mated Animals']

from statistics import mean
from scipy import stats
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

s_bins = 1; fp_fps = 20; post_s = 10; pre_s = 5;

baseline_length = pre_s
event_length = post_s  # FOR ZSCORES

def main():
    global colors, zscore_storage, s_bins, combos, project1;

    rootdir = os.listdir(compile_directory)
    for project_name, groups in zip(project1, combos):
        first = ''
        zscore_storage = make_folder('zscore_summary_byTrial', compile_directory)
        project = make_folder(project_name, zscore_storage)
        zscore_figures = make_folder('Figures', project)
        zscore_bins = make_folder(f"{s_bins}s Bins", project)
        for i in rootdir:
            if i in groups:
                first +=i
                break
        # summary_folder = os.path.join(compile_directory, f"{i}\\Summary") 
        tag_list = compile_directory.split('\\')
        tag_list2 = tag_list[-1].split('_')
        tag = tag_list2[-1]
        # print(tag)
        bin_folder = os.path.join(compile_directory,f'{first}\\Summary\\Behavior_analysis\\during\\zscore_figs\\Averaged per Animal')
        zscore_summaries = os.listdir(bin_folder)

        for zsum in zscore_summaries:
            if '_ALLTRIALS.csv' in zsum:
                for_stats(zsum, bin_folder, groups, first, zscore_bins)
            if '_AVGvalues.csv' in zsum:
                for_figures(zsum, bin_folder, groups, first, zscore_figures)


def for_stats(zsum, bin_folder, groups, first, zscore_bins):
    bin_diary = {}; group_diary = []
    behavior=zsum.split("_ALLTRIALS.csv")
    behavior=behavior[0]
    try:
        zbehavior = pd.read_csv(os.path.join(bin_folder, zsum), header=[0], index_col = None)
        print(behavior)
        bin_log, index =bin_me(zbehavior, s_bins); group_diary.append(first)
        bin_diary.update({f"{first}" : bin_log})
    except pd.errors.EmptyDataError:
        exit()
    for i in groups:
        if i != first:
            bin_folder2 = os.path.join(compile_directory,f'{i}\\Summary\\Behavior_analysis\\during\\zscore_figs\\Averaged per Animal')            
            zscore_summaries2 = os.listdir(bin_folder2)
            if zsum in zscore_summaries2:
                try:
                    zbehavior = pd.read_csv(os.path.join(bin_folder2, zsum), header=[0], index_col = None)
                    bin_log2, _ =bin_me(zbehavior, s_bins); group_diary.append(i)
                    bin_diary.update({f"{i}" : bin_log2})
                except pd.errors.EmptyDataError:
                    continue
    process_bins(bin_diary, group_diary, behavior, index, zscore_bins)
    save_bins(bin_diary, index, zscore_bins, behavior, group_diary)

def for_figures(zsum, bin_folder, groups, first, zscore_figures):
    avg_diary = []; sem_diary = []; group_diary = []
    behavior = zsum.split("_AVGvalues.csv")
    behavior = behavior[0]
    try:
        zbehavior = pd.read_csv(
            os.path.join(bin_folder, zsum), header=[0])
        avg_diary.append(zbehavior["Zscore"]); sem_diary.append(zbehavior["SEM"]); group_diary.append(first)
        fp_times = zbehavior["Time"]
    except pd.errors.EmptyDataError:
        exit()
    for i in groups:
        if i!= first:
            bin_folder2 = os.path.join(compile_directory,f'{i}\\Summary\\Behavior_analysis\\during\\zscore_figs\\Averaged per Animal')            
            zscore_summaries2 = os.listdir(bin_folder2)
            if zsum in zscore_summaries2:
                try:
                    zbehavior2 = pd.read_csv(os.path.join(bin_folder2, zsum), header=[0], index_col = None)
                    avg_diary.append(zbehavior2["Zscore"]); sem_diary.append(zbehavior2["SEM"]); group_diary.append(i)
                except pd.errors.EmptyDataError:
                    continue
    zscores(avg_diary, sem_diary, fp_times,behavior, group_diary, zscore_figures)


def order_groups(groups_):

    while True:
        new_order = input("Would you like to re-order your groups? ")
        if new_order.lower() == 'yes' or new_order.lower() == 'y':
            
            order = list(
                ((input("What's the order? e.g., [1, 0, 2]")).split(',')))
            groups = [groups_[int(i)] for i in order]
            break
        if new_order.lower() == 'no' or new_order == '':
            groups = groups_
            break
    return groups


def save_bins(bin_diary, index, bin_path, behavior, group_diary):
    
    writer = pd.ExcelWriter(f'{bin_path}\{behavior}_bins.xlsx', engine = "xlsxwriter")
    for i in index:
        new_dict = {}
        for g in group_diary:
            new_dict.update({f'{g}': bin_diary[f"{g}"][str(i)]})
        if len(new_dict)!=0:
            
            bins = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in new_dict.items() ]))
            bins.to_excel(writer, sheet_name = str(i))
    writer.close()
def bin_me(zbehavior, s_bins):
    index=-pre_s+s_bins
    index2 = []
    bin_log = []
    range_list  =[i for i in range(0,len(zbehavior), int(fp_fps*s_bins))]
    # print(range_list)
    range_list.append(len(zbehavior))
    bin_log={}
    for i, j in zip(range_list[:-1], range_list[1:]):
        bins = [zbehavior[k][i:j] for k in zbehavior]
        binned = [mean(k) for k in bins]
        bin_log.update({f'{index}': binned})
        index2.append(index)
        index+=s_bins
    return bin_log, index2;

def save_output(tukey_diary, anova_diary, mw_diary, welch, zscore_path2, behavior, index, comp): #https://xlsxwriter.readthedocs.io/example_pandas_conditional.html


    statpath = make_folder("Statistical Summaries", zscore_path2)

    name_lists = []
    df_lists = []
    try:
        if len(welch)!=0:
            writer = pd.ExcelWriter(f'{statpath}\{behavior}_StatsSummary.xlsx', engine = "xlsxwriter")
            Welch = pd.DataFrame(welch, index = index)
            Welch.to_excel(writer, sheet_name = "Unpaired t-test")
            name_lists.append("Unpaired t-test")
            df_lists.append(Welch)
        if len(mw_diary)!=0:
            mw = pd.DataFrame(mw_diary, index = index)
            mw.to_excel(writer, sheet_name = "Mann-Whitney")
            name_lists.append("Mann-Whitney")
            df_lists.append(mw)
        if len(tukey_diary)!=0:
            Tukey = pd.DataFrame(tukey_diary, index = comp)
            Tukey.to_excel(writer, sheet_name = "Tukey")
            name_lists.append("Tukey")
            df_lists.append(Tukey)
        if len(anova_diary)!=0:
            anova = pd.DataFrame(anova_diary, index = None)
            anova.to_excel(writer, sheet_name = "One-Way ANOVA")
            name_lists.append("One-Way ANOVA")
            df_lists.append(anova)
    except ValueError:
        pass
        
    workbook = writer.book
    green = workbook.add_format({'bg_color': 'lime'})
    orange = workbook.add_format({'bg_color': 'orange'})
    yellow = workbook.add_format({'bg_color': 'yellow'})
    blue = workbook.add_format({'bg_color': 'cyan'})

    for i, j in zip(name_lists, range(len(df_lists))):
        t = writer.sheets[str(i)]

        t.conditional_format('A1:P18', {'type': 'cell',
                                                         'criteria': 'between',
                                                         'minimum': 0.0000000000001,
                                                         'maximum': 0.001,
                                                         'format': green,
                                                         })
        t.conditional_format('A1:P18', {'type': 'cell',
                                                         'criteria': 'between',
                                                         'minimum': 0.001000005,
                                                         'maximum': 0.01,
                                                         'format': orange,
                                                         })
        t.conditional_format('A1:P18', {'type': 'cell',
                                                         'criteria': 'between',
                                                         'minimum': 0.010005,
                                                         'maximum': 0.05,
                                                         'format': yellow,
                                                         })
        t.conditional_format('A1:P18', {'type': 'cell',
                                                         'criteria': 'between',
                                                         'minimum': 0.050000005,
                                                         'maximum': 0.099,
                                                         'format': blue,
                                                         })
    writer.close()



def unpaired_ttest(group_diary, bin_diary, index):
    sig_diary = {}
    pairs = []
    for m in range(len(group_diary)) :
        for n in range(len(group_diary)) :
            if m != n:
                pairs.append(tuple((group_diary[m], group_diary[n])))
    pairs = [tuple(sorted(t)) for t in pairs]
    pairs=[*set(pairs)]

    for pair in pairs:
        pvals = []
        for i in index:
            # stat, p_val=stats.ttest_ind(bin_diary[f"{pair[0]}"][i], bin_diary[f"{pair[1]}"][i])
            stat, p_val=stats.ttest_ind(bin_diary[f"{pair[0]}"][str(i)], bin_diary[f"{pair[1]}"][str(i)])
            pvals.append(p_val)
        sig_diary.update({pair: pvals})
    sig_diary_ = pd.DataFrame(sig_diary, index = index)

    return sig_diary_


def manwhit(group_diary, bin_diary, index):
    sig_diary = {}
    pairs = [i for i in zip(group_diary[:-1], group_diary[1:])]
    pairs = [tuple(sorted(t)) for t in pairs]
    pairs=[*set(pairs)]

    for pair in pairs:
        pvals = []
        for i in index:
            stat, p_val=stats.mannwhitneyu(bin_diary[f"{pair[0]}"][str(i)], bin_diary[f"{pair[1]}"][str(i)])
            
            pvals.append(p_val)
        sig_diary.update({pair: pvals})
    sig_diary_ = pd.DataFrame(sig_diary, index = index)

    return sig_diary_

def tukey_test(bin_diary, group_diary, m):
    pairs = []
    pval = []
    try:
        if len(group_diary)>2:
            if len(group_diary)==4:
                res = stats.tukey_hsd(bin_diary[f"{group_diary[0]}"][str(m)], bin_diary[f"{group_diary[1]}"][str(m)], bin_diary[f"{group_diary[2]}"][str(m)], bin_diary[f"{group_diary[3]}"][str(m)])
            if len(group_diary)==3:
                res = stats.tukey_hsd(bin_diary[f"{group_diary[0]}"][str(m)], bin_diary[f"{group_diary[1]}"][str(m)], bin_diary[f"{group_diary[2]}"][str(m)])
            if len(group_diary)==5:
                res = stats.tukey_hsd(bin_diary[f"{group_diary[0]}"][str(m)], bin_diary[f"{group_diary[1]}"][str(m)], bin_diary[f"{group_diary[2]}"][str(m)], bin_diary[f"{group_diary[3]}"][str(m)], bin_diary[f"{group_diary[4]}"][str(m)])
            #DUNNETT REQUIRES A CONTROL
            # stat, p_val = stats.dunnett(bin_diary[f"{group_diary[0]}"][str(i)], bin_diary[f"{group_diary[1]}"][str(i)], bin_diary[f"{group_diary[2]}"][str(i)], bin_diary[f"{group_diary[3]}"][str(i)], control = controls)
            for ((i, j), l) in np.ndenumerate(res.pvalue):
                if i!=j:
                    pval.append(l)
                    g = str(group_diary[i])
                    h = str(group_diary[j])
                    pairs.append((g, h))
            # pairs = [tuple(sorted(t)) for t in pairs]
    except IndexError:
        pass
    except ValueError:
        pass
    # print(len(pairs), pairs)

    return pval, pairs


def process_bins(bin_diary, group_diary, behavior, index, zscore_path2):

#TUKEY
    tukey_diary = {}
    anova_vals = []
    anova_diary = {}
    for m in index:
        try:
            pval, comp = tukey_test(bin_diary, group_diary, m)
            # print(behavior, bin_diary[f"{group_diary[0]}"][str(m)])
            try:
                if len(group_diary)==4:
                    stat, anova = stats.f_oneway(bin_diary[f"{group_diary[0]}"][str(m)], bin_diary[f"{group_diary[1]}"][str(m)], bin_diary[f"{group_diary[2]}"][str(m)], bin_diary[f"{group_diary[3]}"][str(m)])
                elif len(group_diary)==3:
                    stat, anova = stats.f_oneway(bin_diary[f"{group_diary[0]}"][str(m)], bin_diary[f"{group_diary[1]}"][str(m)], bin_diary[f"{group_diary[2]}"][str(m)])
                else:
                    continue
            except IndexError or ValueError:
                pass
            anova_vals.append(anova)   
            for (a,b), z in zip(comp, pval):
                for (w,x), y in zip(comp, pval):
                    if (a,b)==(x,w):

                        comp.remove((a,b))
                        pval.remove(z)
                        break
        except ValueError:
            pass
        tukey_diary.update({f'{m}': pval})
    anova_diary.update({"Index": index})
    anova_diary.update({"P-value": anova_vals})


    mw_diary = manwhit(group_diary, bin_diary, index)
    ttest =unpaired_ttest(group_diary, bin_diary, index)
    save_output(tukey_diary, anova_diary, mw_diary, ttest, zscore_path2, behavior, index, comp)

def zscores(avg_diary, sem_diary, fp_times, j, group_diary, zscore_path3):
    print(j)

    color_me=[colors[f"{d}"] for d in group_diary]

    zfigures=plt.figure()
    mpl.rcParams['axes.linewidth'] = 1.3

    plt.title(f"{j}", y = 1.05)
    for i in range(len(avg_diary)):
        plt.plot(fp_times, avg_diary[i], label=group_diary[i], color=color_me[i], alpha=1)
        plt.fill_between(fp_times, avg_diary[i]-sem_diary[i], avg_diary[i] +
                            sem_diary[i], alpha=0.1, color=color_me[i])
    plt.xticks(np.arange(-baseline_length, event_length+1), size = 15)
    plt.yticks(size = 15)
    # plt.axvline(x=0, color='k', alpha=0.5)
    plt.legend(loc='upper left', fontsize = 12)
    plt.xlabel("Time (s)", size = 15)
    plt.ylabel("Z-score", size = 15)
    plt.xlim(-5,10)

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    # plt.ylim(-5, 15)
    
    # plt.show()

    zfigures.savefig(f'{zscore_path3}\{j}.tif')
    zfigures.savefig(f'{zscore_path3}\{j}.svg')


    plt.close()

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





if __name__ == "__main__":  # avoids running main if this file is imported
    main()