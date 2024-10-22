# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
# https://stacbehav_filesoverflow.com/questions/11190735/python-matplotlib-superimpose-scatter-plots

# project_home="F:\Photometry-Fall2022\Pup Block\Virgin\Trial 1-PE Males\Ball"
# PROCEED WITH CAUTION CHANGING ANYTHING AFTER THIS LINE
# MAKE SURE ALL PACKAGES ARE INSTALLED

remove_traces=[] #NOTE only input their id, separate by commas

colors="red"

# project_home = "E:\Photometry-Fall2022\Final Analysis\REANALYSIS_zscore_nolimit\Virgin Females\Pup"

from combined_photocode_scaled_trace_05052024 import make_folder, project_home, lolo, fp_fps, post_s, pre_s
import pandas as pd
import numpy as np
import os
import csv
from scipy import stats
import matplotlib.pyplot as plt
from statistics import mean
# import seaborn as sns


baseline_length=pre_s
event_length=post_s #FOR ZSCORES

root = os.path.join(project_home, "Behaviors")
analyses = os.listdir(root)

summary_path=make_folder("Summary", project_home)

def main():

    for analysis in analyses:
        if analysis=="Behavior_analysis":
            analyze=os.path.join(root, analysis)
            states=os.listdir(analyze)
            for state in states: 
                if state=='during' or state=='controls':
                    print(state.upper())
                    entrypoint = os.path.join(analyze, state)
                    behaviors = os.listdir(entrypoint)
                    # CREATES SUB-FOLDER FOR EACH STATE IN SUMMARY FOLDER
                    analysis_folder=make_folder(analysis, summary_path)
                    state_save=make_folder(state, analysis_folder)
                
                    percentage_diary = {}; perc_avg={}
                    zscore_avg = {}; auc_zscore={}
                    mag_auc={}; mag_auc_avg={};
                    auc_diary = {}; auc_avg={}
                    raw_diary_i = {}; raw_diary_avg = {}
                    amplitude = {}; amplitude_avg={}
                    deltaF = {}; delf_avg={}; all_zscore={}
                    probability_peak = {}; probability_behav = {}; zscore_ind={}
                    behavior_list = []
                    # if state!="before":

                    # ACCESS FOLDERS OF BEHAVIORS
                    for i in behaviors:
                        print(i)
                        _ = os.path.join(entrypoint, i)
                        behav_files = os.listdir(_)
                        behavior_list.append(i)

                        # VARIABLES TO SAVE FOR EACH BEHAVIOR; EMPTY WHEN NEXT BEHAVIOR FOLDER OPENS
                        event_means = []  # Decide if I wanna append or extend
                        baseline_means = []
                        entrance_event = []; entrance_baseline = []
                        tag_id1 = []
                        meansof_baseline = []; meansof_event = []
                        peri_event_matrix = []; peri_baseline_matrix = []
                        entrance_base_matrix = []; entrance_event_matrix = []
                        prominences = []; prom_avg =[]
                        bases = []; place_averages=[]
                        delF = []; del_avgg=[]
                        peak_prob = []; behav_prob = []
                        avg_matrix_base=[]; avg_matrix_event=[]
                        
                        # ACCESS SAVED VARIABLES IN EACH BEHAVIOR FOLDER
                        for k in behav_files:
                            # EACH INDIVIDUAL EVENT MEANS
                            if 'event_means' in k:
                                l = os.path.join(_, k)
                                q = k.split('_')
                                r = q[-1].split('.')
                                tag = r[0]  # make sure the tag include _#.CSV and not #.CSV
                                # print(tag)
                                # if str(tag) in remove_traces:
                                #     print(tag)
                                #     continue
                                event = list(csv.reader(open(l)))
                                event_means.extend(event)

                                behav_files2=behav_files
                                
                                # FIND MATCHING BASELINE FILE
                                for j in behav_files2:
                                    q = j.split('_')
                                    r = q[-1].split('.')
                                    tag2=r[0]
                                    if 'baseline_means' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        base = list(csv.reader(open(l)))
                                        baseline_means.extend(base)
                                        pp = [str(tag) for i in range(len(base))]
                                        tag_id1.extend(pp)

                                        if i == lolo and state!='before':
                                            entrance_event.extend(event[0])
                                            behavior_list.append('Entrance')
                                            entrance_baseline.extend(base[0])

                                    if 'peri_event_matrix' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        event_matrix = list(csv.reader(open(l)))
                                        avg_matrix_event.append(np.average(np.asarray(event_matrix, dtype=float), axis=0))

                                        peri_event_matrix.extend(event_matrix)

                                        if i == lolo and state!='before':
                                            entrance_event_matrix.append(event_matrix[0])

                                    if 'peri_baseline_matrix' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        baseline_matrix = list(csv.reader(open(l)))
                                        avg_matrix_base.append(np.average(np.asarray(baseline_matrix, dtype=float), axis=0))
                                        peri_baseline_matrix.extend(baseline_matrix)
                                        if i == lolo and state!='before':
                                            entrance_base_matrix.append(baseline_matrix[0])

                                    if 'place_value' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        place_value = list(csv.reader(open(l)))
                                        place_averages.append(mean(buoyant(flatten(place_value))))
                                        bases.extend(place_value)
                                    if 'peak_prom' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        proms = list(csv.reader(open(l)))
                                        prom_avg.append(mean(buoyant(flatten(proms))))
                                        prominences.extend(proms)
                                    if 'amplitude' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        amp = list(csv.reader(open(l)))
                                        del_avgg.append(mean(buoyant(flatten(amp))))
                                        delF.extend(amp)
                                    if 'probability' in j and tag==tag2:
                                        l = os.path.join(_, j)
                                        peep = list(csv.reader(open(l)))
                                        peak_prob.extend(peep[0])
                                        behav_prob.extend(peep[1])

                                    if i != lolo:
                                        # AVERAGED EVENT MEAN PER ANIMAL
                                        if 'means_to_compare' in j and tag==tag2:
                                            l = os.path.join(_, j)
                                            csvfile = list(csv.reader(open(l)))
                                            meansof_event.extend(csvfile[1])
                                            # need to drop the index and details
                                            meansof_baseline.extend(csvfile[0])

                        if len(baseline_means)!=len(event_means):
                            exit("Length mismatch. Check code; Individual values will be affected")
                        if i == lolo:# and state!='before':
        #percent
                            _, entrance_percent = percent_find(
                                i, entrance_event, entrance_baseline)
                            u = "entrance"
                            percentage_diary[f"ANIMAL ID_Entrance"] = tag_id1
                            percentage_diary["Entrance baseline"] = [
                                100 for i in range(len(entrance_percent))]
                            percentage_diary["Entrance"] = entrance_percent

                            perc_avg["ANIMAL ID"] = list(dict.fromkeys(tag_id1))
                            perc_avg["Entrance baseline"] = [
                                100 for i in range(len(entrance_percent))]
                            perc_avg["Entrance"] = entrance_percent

                        if len(peri_baseline_matrix) == len(peri_event_matrix): # and len(peri_baseline_matrix)!=0 and i!=lolo:
        #INDV Zscore 
                #INDV Zscore 
                            try:
                                u="Single Event"
                                base_avg, event_avg, event_stacks, base_stacks, diary2, path=zscores(peri_baseline_matrix, peri_event_matrix, i, state_save, u, tag_id1)
                                all_zscore[f"ANIMAL ID_{i}"]=tag_id1
                                # full_zplot(diary, i, path)
                                zscore_ind[f"ANIMAL ID_{i}"]=tag_id1
                                # zscore_avg[f"{i} baseline"] = base_avg
                                zscore_ind[f"{i}"] = [e-b for e,b in zip(event_avg,base_avg)]
                            except UnboundLocalError or NotADirectoryError:
                                continue
#AVG Zscore        
                            try:
                                u="Averaged per Animal"
                                base_avg, event_avg, event_stacks, base_stacks, diary2, path=zscores(avg_matrix_base, avg_matrix_event, i, state_save, u, list(dict.fromkeys(tag_id1)))
                                all_zscore[f"ANIMAL ID_{i}"]=tag_id1
                                # full_zplot(diary, i, path)
                                zscore_avg[f"ANIMAL ID_{i}"]=list(dict.fromkeys(tag_id1))
                                # zscore_avg[f"{i} baseline"] = base_avg
                                zscore_avg[f"{i}"] = [e-b for e,b in zip(event_avg,base_avg)]
                            except UnboundLocalError:
                                continue
        #AUC ZSCORE
                            auc_event, auc_base = area_under_the_curve(
                                event_stacks, base_stacks)
                            auc_zscore[f"ANIMAL ID_{i}"] = tag_id1
                            # auc_zscore[f"{i} baseline"] = auc_base
                            auc_zscore[f"{i}"] = [e-b for e,b in zip(auc_event,auc_base)]

        #Regular AUC
                            auc_event, auc_base = area_under_the_curve(
                                peri_event_matrix, peri_baseline_matrix)
                            auc_diary[f"ANIMAL ID_{i}"] = tag_id1
                            # auc_diary[f"{i} baseline"] = auc_base
                            auc_diary[f"{i}"] = [e-b for e,b in zip(auc_event,auc_base)]    
        #AVG AUC
                            avg_event_auc, avg_base_auc = area_under_the_curve(avg_matrix_event, avg_matrix_base)
                            auc_avg[f"ID_{i}"] = list(dict.fromkeys(tag_id1))
                            # auc_avg[f"{i} baseline"]=avg_base_auc
                            auc_avg[f"{i}"]= [e-b for e,b in zip(avg_event_auc, avg_base_auc)]
        #AUC MAG
                            # factor=[100/m for m in auc_base]
                            auc_event=[((10/m)*n-10) for m, n in zip(auc_base,auc_event)]
                            mag_auc[f"ANIMAL ID_{i}"] = tag_id1
                            # mag_auc[f"{i} baseline"] = auc_base
                            mag_auc[f"{i}"] = auc_event     
        #AVG AUC MAG
                            # factor=[100/m for m in avg_base_auc]
                            # avg_matrix_base=[100 for _ in range(len(factor))]
                            mag_auc_avg[f"ANIMAL ID_{i}"] = list(dict.fromkeys(tag_id1))
                            avg_event_auc=[((10/m)*n-10) for m, n in zip(avg_base_auc,avg_event_auc)]

                            # mag_auc_avg[f"{i} baseline"]=avg_matrix_base
                            mag_auc_avg[f"{i}"]= avg_event_auc

        #Raw Values
                        # RAW VALUES- INDIVIDUAL

                        # baseline_means=[float(i) for i in baseline_means]
                        if len(baseline_means)!=0 and i!=lolo:
                            raw_diary_i[f"ANIMAL ID_{i}"] = tag_id1
                            # raw_diary_i[f"{i} baseline"] = flatten(baseline_means)
                            raw_diary_i[f"{i}"] = [float(e)-float(b) for e,b in zip(flatten(event_means), flatten(baseline_means))]
                            raw_diary_i['        '] = []

                            # RAW VALUES AVERAGE
                            average_base_matrix_list = np.asarray([mean(i[(pre_s*fp_fps):]) for i in avg_matrix_base])
                            average_event_matrix_list = np.asarray([mean(i) for i in avg_matrix_event])

                            # print("final length: ", len(average_event_matrix_list))

                            raw_diary_avg[f"ANIMAL ID_{i}"] = list(dict.fromkeys(tag_id1))
                            # raw_diary_avg[f"{i} baseline"] = average_base_matrix_list
                            raw_diary_avg[f"{i}"] = [e-b for e,b in zip(average_event_matrix_list,average_base_matrix_list)]
                            raw_diary_avg['           '] = []

        #Percent Change
                        # Percent_Baseline for each behavior
                        if len(event_means)!=0 and i!=lolo:
                            # print(event_means)
                            percent_baseline, percent_event = percent_find(
                                i, event_means, baseline_means)
                            _, perc_event_avg=percent_find(i, average_event_matrix_list, average_base_matrix_list)
                            # print(average_base_matrix_list)
                            if len(percent_event)!=0 and i!=lolo:
                                percentage_diary[f"ANIMAL ID_{i}"] = list(dict.fromkeys(tag_id1))
                                # percentage_diary[f"{i} baseline"] = [
                                #     100 for i in range(len(percent_event))]
                                percentage_diary[f"{i}"] = percent_event

                                perc_avg[f"ANIMAL ID_{i}"] = list(dict.fromkeys(tag_id1))
                                perc_avg[f"{i}"] = perc_event_avg


                        average_sums = make_folder("Averages", state_save)

                        indv_vals = make_folder("Individual Values", state_save)
                    o="Percentage_Change"
                    # if analysis=="Peak_analysis":
                    #     save_dict(percentage_diary, indv_vals, "Percentage_Change")    
                    #     save_dict(auc_diary, indv_vals, "AUC") 
                    #     save_dict(raw_diary_i, indv_vals, "Individual Raw Means")
                    #     save_dict(amplitude,indv_vals, "Amplitude")
                    #     save_dict(raw_diary_avg, average_sums, "Average Raw Means")
                    #     save_dict(deltaF,indv_vals, "DeltaF")
                    #     save_dict(probability_peak, average_sums, "Peak Probability")
                    #     save_dict(probability_behav, average_sums, "Behavior Probability")
                    #     save_dict(perc_avg, average_sums, "Average Percent Change")
                    #     save_dict(auc_avg, average_sums, "Average AUC")
                    #     save_dict(delf_avg, average_sums, "Average DeltaF")
                    #     save_dict(amplitude_avg, average_sums, "Average Amplitude")
                    #     save_dict(zscore_avg, average_sums, "Average Z-score")
                    #     # save_dict(auc_zscore, average_sums, "Average AUC ZSCORE")
                    #     save_dict(mag_auc, indv_vals, "Magnitude of AUC")
                    #     save_dict(mag_auc_avg, average_sums, "Average Magnitude AUC")
                    if analysis=="Behavior_analysis":

                        save_dict(percentage_diary, indv_vals, "Percentage_Change")    
                        save_dict(auc_diary, indv_vals, "AUC") 
                        save_dict(raw_diary_i, indv_vals, "Individual Raw Means")
                        save_dict(raw_diary_avg, average_sums, "Average Raw Means")
                        save_dict(perc_avg, average_sums, "Average Percent Change")
                        save_dict(auc_avg, average_sums, "Average AUC")
                        save_dict(zscore_ind, indv_vals, "Individual Z-score")
                        save_dict(zscore_avg, average_sums, "Average Z-score")
                        save_dict(auc_zscore, average_sums, "Average AUC ZSCORE")
                        save_dict(mag_auc, indv_vals, "Magnitude of AUC")
                        save_dict(mag_auc_avg, average_sums, "Average Magnitude AUC")


        print("Merger Complete")

def buoyant(lst):
    x= [float(i) for i in lst]
    return x

def save_dict(x,w, o):
    pc = pd.DataFrame(dict([(k, pd.Series(v))
                        for k, v in x.items()]))
    pc.to_csv(
        f"{w}/{o}_summary.csv", header=True)
    
def area_under_the_curve(peri_event_matrix, peri_baseline_matrix):
    event = []
    base = []
    peri_baseline = np.asarray(peri_baseline_matrix, dtype=float)
    peri_event = np.asarray(peri_event_matrix, dtype=float)

    for i, j in zip(peri_event, peri_baseline):
        i.astype(float)
        j.astype(float)
        min_base = j.min()
        auc_event = np.abs(np.trapz(i))
        auc_base = np.abs(np.trapz(j))
        event.append(auc_event)
        base.append(auc_base)

    return event, base


def flatten(x):
    list = [i for ii in x for i in ii]
    return list


def percent_find(i, event_means, baseline_means):

    base = []
    event = []

    event_means = np.array(event_means)
    baseline_means = np.array(baseline_means)

    for x, y in zip(event_means, baseline_means):
        x = float(x)
        y = float(y)
        # print(x, y, i)
        a = y/y*100
        b = ((x-y)/y)*100  # event/baseline*100
        base.append(a)
        event.append(b)
    base = np.array(base)
    event = np.array(event)
    return base, event


def zscores(peri_baseline_matrix, peri_event_matrix, i, state_save, u, tag_id):

    zscore_diary = []
    fp_times=[]
    zscore_diary2=[]

    peri_base = np.asarray(peri_baseline_matrix, dtype=float)
    peri_event = np.asarray(peri_event_matrix, dtype=float)

    ##Individual files

    #average for each animal
    try:
        for o, z in zip(peri_base, peri_event):
            # a=np.average(o[int(len(o)/2):])
            # b=np.std(o[int(len(o)/2):])
            #5 SECOND BASELINE
            a=np.average(o[int(-baseline_length*fp_fps):])
            b=np.std(o[int(-baseline_length*fp_fps):])
            zscore=[(c-a)/b for c in o[int(-baseline_length*fp_fps):]]
            zscore.extend([(c-a)/b for c in z])
            zscore_diary.append(zscore)

            #10 SECOND BASELINE
            a1=np.average(o[int(-event_length*fp_fps):])
            b1=np.std(o[int(-event_length*fp_fps):])
            zscore2=[(c-a1)/b1 for c in o[int(-event_length*fp_fps):]]
            zscore2.extend([(c-a1)/b1 for c in z])
            zscore_diary2.append(zscore2)


        p = np.arange(int(-baseline_length*fp_fps), 0); low = p/fp_fps
        q = np.arange(int(event_length*fp_fps)); high = q/fp_fps

        fp_times.extend(low); fp_times.extend(high)

        base_avg=[np.mean(f[:int(baseline_length*fp_fps)]) for f in zscore_diary]
        event_avg=[np.mean(f[int(baseline_length*fp_fps):]) for f in zscore_diary]
        # stat, p_val=stats.ttest_1samp(event_avg, popmean=0)
        stat, p_val=stats.ttest_rel(base_avg, event_avg,)
        diary = pd.DataFrame(zscore_diary).T  # Modified line        # sns.heatmap(zscore_diary, xticklabels=fp_times)
        diary.columns = tag_id 

        averages=np.average(zscore_diary, axis=0)
        sem=stats.sem(zscore_diary, axis=0)
        # std=stats.std(zscore_diary, axis=0)

        save_it=pd.DataFrame({"Zscore": np.asarray(averages), "SEM": np.asarray(sem), "Time": np.asarray(fp_times)})

        # t=(base_avg-event_avg)/std*np.sqrt((1/len(peri_base))+(1/len(peri_event)))

        zscore_fig = plt.figure()
        plt.plot(fp_times, averages, alpha=1, color=colors)
        plt.ylabel('Z-Score'); plt.xlabel('Time (s)')
        plt.xticks(np.arange(-baseline_length, event_length))
        plt.title(f"{i}_p={p_val:.3f}")  # fix y labels and tighten the graph
        plt.fill_between(fp_times, averages-sem, averages +
                            sem, alpha=0.1, color=colors)
        plt.axvline(x=0, color='k')
        # plt.ylim((-6, 20))

        # plt.show()
        
        zscore_path = make_folder("zscore_figs", state_save)
        zscore_path = make_folder(u, zscore_path)

        zscore_fig.savefig(f'{zscore_path}\{i}_zscore_figs.tif')
        zscore_fig.savefig(f'{zscore_path}\{i}_zscore_figs.svg')
        # plt.close()

        # labels=[round(p) for p in fp_times]
        # df=pd.DataFrame(zscore_diary, columns=labels)
        # heat=plt.figure()
        # sns.heatmap(df, cmap='coolwarm', center=1, robust=True, xticklabels=(pre_s+post_s+1))
        # plt.axvline(x=(pre_s*fp_fps), color='k')
        # plt.title(f"{i}")
        # heat.savefig(f'{zscore_path}\{i}_HEATMAP.tif')

# def full_zplot(zscore_diary, i, zscore_path):
#     # m=int((len(zscore_diary[0])/2)/fp_fps)
#     # labels=[round(p) for p in np.arange(-m,m,0.05)]
#     labels=[round(p) for p in np.arange(-pre_s, post_s,0.05)]
#     df=pd.DataFrame(zscore_diary, columns=labels)
#     heat=plt.figure()
#     sns.heatmap(df, cmap='coolwarm', center=1, robust=True, xticklabels=((pre_s+post_s)))
#     plt.title(f"{i}")
#     heat.savefig(f'{zscore_path}\{i}_HEATMAP.tif')

        save_it.to_csv(f'{zscore_path}\{i}_AVGvalues.csv', index=False, header=True)
        diary.to_csv(f'{zscore_path}\{i}_ALLTRIALS.csv', index=False, header=True)
        plt.close()

        event_stacks=[zscore[int(event_length*fp_fps):] for zscore in zscore_diary2]
        base_stacks=[zscore[0:int(event_length*fp_fps)] for zscore in zscore_diary2]

    except IndexError:
        print("error")
        exit
    except ValueError:
        print("error")
        exit
    # zscore_diary = np.asarray(zscore_diary, dtype=float)        

    return base_avg, event_avg, event_stacks, base_stacks, zscore_diary, zscore_path;



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

def plot_dF_trace(peri_baseline_matrix, peri_event_matrix, i, state_save, u):

    dF_diary = []
    fp_times=[]
    peri_base = np.asarray(peri_baseline_matrix, dtype=float)
    peri_event = np.asarray(peri_event_matrix, dtype=float)

    #average for each animal
    try:
        for o, z in zip(peri_base, peri_event):
            baseline_end_index = int(baseline_length*fp_fps)
            dF = [c for c in o[-baseline_end_index:]]  # Z-score calculation                    
            dF.extend([c for c in z])
            dF_diary.append(dF)

        p = np.arange(int(-baseline_length*fp_fps), 0); low = p/fp_fps
        q = np.arange(int(event_length*fp_fps)); high = q/fp_fps

        fp_times.extend(low); fp_times.extend(high)

        base_avg=[np.mean(f[:int(baseline_length*fp_fps)]) for f in dF_diary]
        event_avg=[np.mean(f[int(baseline_length*fp_fps):]) for f in dF_diary]
        # stat, p_val=stats.ttest_1samp(event_avg, popmean=0)
        stat, p_val=stats.ttest_rel(base_avg, event_avg,)
        diary = pd.DataFrame(dF_diary).T
        # sns.heatmap(zscore_diary, xticklabels=fp_times)
        # plt.show()

        averages=np.average(dF_diary, axis=0)
        sem=stats.sem(dF_diary, axis=0)
        # std=stats.std(zscore_diary, axis=0)
        save_it=pd.DataFrame({"Zscore": np.asarray(averages), "SEM": np.asarray(sem), "Time": np.asarray(fp_times)})

        # t=(base_avg-event_avg)/std*np.sqrt((1/len(peri_base))+(1/len(peri_event)))
        df_trace = plt.figure()
        plt.plot(fp_times, averages, alpha=1, color=colors)
        plt.ylabel('zF'); plt.xlabel('Time (s)')
        plt.xticks(np.arange(-baseline_length, event_length))
        plt.title(f"{i}_p={p_val:.3f}")  # fix y labels and tighten the graph
        plt.fill_between(fp_times, averages-sem, averages +
                            sem, alpha=0.1, color=colors)
        plt.axvline(x=0, color='k')
        # plt.ylim((-6, 20))

        # plt.show()
        
        dF_path = make_folder("trace_ROI", state_save)
        dF_path = make_folder(u, dF_path)

        df_trace.savefig(f'{dF_path}\{i}_ROI_figs.tif')
        df_trace.savefig(f'{dF_path}\{i}_ROI_figs.svg')
        # plt.close()

        # labels=[round(p) for p in fp_times]
        # df=pd.DataFrame(zscore_diary, columns=labels)
        # heat=plt.figure()
        # sns.heatmap(df, cmap='coolwarm', center=1, robust=True, xticklabels=(pre_s+post_s+1))
        # plt.axvline(x=(pre_s*fp_fps), color='k')
        # plt.title(f"{i}")
        # heat.savefig(f'{zscore_path}\{i}_HEATMAP.tif')

        save_it.to_csv(f'{dF_path}\{i}_AVGvalues.csv', index=False, header=True)
        diary.to_csv(f'{dF_path}\{i}_ALLTRIALS.csv', index=False, header=True)
        plt.close()


    except IndexError:
        print("error")
        exit
    except ValueError:
        print("error")
        exit

if __name__ == "__main__":  # avoids running main if this file is imported
    main()



