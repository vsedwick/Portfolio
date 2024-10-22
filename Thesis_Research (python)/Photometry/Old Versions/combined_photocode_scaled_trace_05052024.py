# AUTHOR: VICTORIA SEDWICK
# ADAPTED FROM NEUROPHOTOMETRICS AND ILARIA CARTA (MATLAB)

# https://mat#plotlib.org/stable/gallery/sub#plots_axes_and_figures/sub#plots_demo.html
# https://stackoverflow.com/questions/11190735/python-mat#plotlib-superimpose-scatter-#plots

# <Where are your files located?>

project_home = r"E:\VS-Jade_Mated_odor_photometry\scaled_trace\Mated Females"
# project_home=input("Project folder?: ")
# NOTE <Load any ethovisiol
# <Frames per second for photo and ethovision data>n file that has behavior names listed>
example_behav_file = r"E:\Photometry-Fall2022\Final Analysis\REANALYSIS_zscore_nolimit\V. Males Inf\Pups\45\Raw data-Pup Exposure 2nd set-Trial    38 (2).xlsx"

behav_fps = 30
fp_fps = 20
# <How many seconds do you want to analyze for peri-event?>
pre_s = 5; post_s = 10; offset = 0
#if full_trace = yes
crop_end=2; crop_front=2
duration_s=0.1; behavior_row=35
#YES OR NO
35
full_trace='no'  #if no, the program will crop around the start and end times
use_zones='yes'
use_groups='yes'
controls='yes'  

#only if full_trace=no
time_from_start=120; time_from_end=120
    
stop_at=0
# <Which behaviors would be categorized as pup interaction?>
# cotton_behaviors = ['qtip']
# stick_behaviors= ['stick']
stim_interaction_all = []
# non_stim_interaction = ['Rearing', 'Digging']
# min_duration = ["Digging", "Rearing", "Grooming", 'Retrieval']
# aggression = ["Attack", "Aggressive groom"]   
stim_contact = ['Grooming', 'Nudge']
stim_interaction= ['Grooming', 'Sniff', 'Nudge'] 
non_stim_interaction = ['Rearing', 'Digging']
# stim_contact = ['Grooming', 'Nudge']
# stim_interaction= ['Grooming', 'Sniff', "Nudge"] 
# non_stim_interaction = ['Rearing', 'Digging']
min_duration = ["Digging", "Rearing", "Grooming"]


objects=['qtip', 'stick']
check_in_zone=['Sniff', 'Grooming', 'Digging', 'Rearing']
sniff_group=['sniff', 'Grooming']


zone = str('In nest')

# Which roi/fiber signal do you want to #plot (0-2)?
z = 0
# name of start behavior
lolo = str('start')
first = str('Approach')
point_events=['Approach', 'start', 'bite']

Groups = ({'Stim Interaction':stim_interaction,
           'Stim Contact': stim_contact, 
           'No Stim Interaction': non_stim_interaction
        #    'Aggressive Behaviors': aggression
        })

# Groups = ({ mke
#            'Cotton tip': cotton_behaviors,
#            "Stick int": stick_behaviors,
#            'No Stim Interaction': non_stim_interaction,
#            'Stim Interaction': stim_interaction})

# NOTE PROCEED WITH CAUTION CHANGING ANYTHING AFTER THIS LINE
# MAKE SURE ALL PACKAGES ARE INSTALLED WITH 'pip install [package]'
import matplotlib.colors
import tkinter
import os
import csv
from math import floor
from numpy import mean
from scipy import stats, signal
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import json


def main():
    global behav_fps, behavior_row, fp_fps, example_behav_file, z, pre_s, post_s, lolo, controls, Groups, use_zones, use_groups, zone, crop_front, duration_s, min_duration, point_events, first, contact_behaviors, stop_at, complex_behaviors;

    for _ in str(project_home):
        if _ == '\\':
            project_home.replace(_, '\\\\')
    print(project_home)

    # LOAD EXAMPLE BEHAVIOR FILE FOR LIST OF AVAILABLE BEHAVIORS
    score_behavior = open_gui(example_behav_file, behavior_row)
    report=[]

# NAVIGATE THROUGH DATA FOLDERS
    root = os.listdir(project_home)
    for animal in root:
        ext=[".avi",".csv",".mp4", ".pzfx", ".jpg", ".txt", ".prism", ".xlsx", ".doc"]
        if animal == "Behaviors_PA" or animal == "Behaviors" or animal == "Summary" or animal == "Summary_PA" or animal == "Archive" or animal == "Videos" or animal.endswith('.CSV') or animal.lower().endswith(tuple(ext)):
            continue
        else:
            project_id = animal
            animal_path = os.path.join(project_home, animal)
            file_to_load = os.listdir(animal_path)

            # LOAD FILESpre_
            print(f'Analyzing ID: {animal}...')
            print(f"{animal}. Opening files")
            try:
                fp_raw, behav_raw, report = open_files(file_to_load, animal_path, animal, report)
            except UnboundLocalError:
                continue

            fp_time, _470, _410 = fp_split(fp_raw, z, fp_fps) 
            behavior, behav_frames, behav_time = behav_split(behav_raw, behav_fps)

            ##PAD PHOTOFILES
            # if round(len(_470)/behav_time.max())!= fp_fps:
            #     if int(len(_470)/behav_time.max()) < fp_fps-4 or int(len(_470)/behav_time.max())>fp_fps:
            #         report.append(animal)
            #         continue
            #     else:
            print(f"Correcting frames per second. Current fps: {int(len(_470)/behav_time.max())}")
            _470, _410, fp_time, fp_frames = correct_photo_fps(_470, _410, behav_time.max(), fp_time.max())
            print(f"New frames per second: {int(len(_470)/fp_time.max())}")

            print(f"Time difference between video and trace before cropping: {behav_time.max()}, {fp_time.max()}")

            #ALIGN BEHAVIOR AND PHOTOMETRY

            #crops the front
            behav_time, behav_raw, cut=takeoff(behav_time, fp_time, behavior, behav_fps, behav_raw)  #wanna crop end time (100?)
            start_times=np.array(behav_raw[lolo])
            place=[int(i) for i in range(len(start_times)) if start_times[i]==1]
            try:
                if controls=='yes': 
                    if len(place)<4:
                        start=place[0]; end=place[-1];
                    else:
                        start=place[2]; end=place[3];
                
                elif controls=='no': start=place[0]; end=place[1]
            except IndexError:
                print("The start parameter is not indicated for this subject")
                report.append(animal)
                continue
            
            fp_start=floor(start/behav_fps*fp_fps); fp_end=floor(end/behav_fps*fp_fps)
            if full_trace.lower()=='no':
                _470,_410, behav_raw, behav_time, fp_time, fp_frames, behav_frames=crop_trace(_470, _410, behav_raw, behavior, place[0], place[-1])
            else:           
                _470,_410, behav_raw, behav_time, fp_time, fp_frames, behav_frames=cropendtime(_470, _410, behav_raw, behavior, crop_end) #cuts the last minute

            print(f"Time difference between video and trace after cropping: {behav_time.max()}, {fp_time.max()}")

            print(f'{animal}. Normalizing photometry trace')
            smooth_trace, fp_frames, heights, l_base=photo_info(_470, _410, fp_frames, project_id)       
            if type(smooth_trace)==str and smooth_trace.lower() == 's':
                report.append(animal)
                continue     
            #MAKE AS BOOLEAN ARRAY           
            if use_zones=='yes': 
                behav_scores, score_behaviors=zones_in_use(check_in_zone, behav_raw, score_behavior, zone)
                print(score_behaviors)
            if use_zones=='no':
                score_behaviors=score_behavior
                behav_scores={f"{i}": np.array(behav_raw[i], dtype=bool) for i in score_behaviors}
                # print(score_behaviors)
            behavior_paths=make_folder("Behaviors", project_home)

            #START TIMES & CONTROLS OR NO CONTROLS
            start_times=np.array(behav_raw[lolo])
            place=[int(i) for i in range(len(start_times)) if start_times[i]==1]
            start_placements = [1]
            end_placements = [place[0]]
            stages = ['before', 'during']

            if len(place)<4 and controls=='yes':
                stages[0] = 'controls'
                print('len is <4')
                start_placements.append(place[0]); end_placements.append(place[-1])
            elif controls=='yes': 
                stages.append('controls')
                #DURING
                start_placements.append(place[2]); end_placements.append(place[-1])
                #Controls
                start_placements.append(place[0]); end_placements.append(place[1])
            else:
                start_placements.append(place[0]); end_placements.append(place[-1])

            peak_f=make_folder("Peak_analysis", behavior_paths); behav_f=make_folder("Behavior_analysis", behavior_paths);
            total_fpstart = int(place[0]/behav_fps*fp_fps); total_fpend= int(place[-1]/behav_fps*fp_fps)
            
            for start,end,stage in zip(start_placements,end_placements,stages):            
                fp_start=int(start/behav_fps*fp_fps); fp_end=int(end/behav_fps*fp_fps)
                peak_folder=make_folder(stage, peak_f); behav_folder=make_folder(stage, behav_f)
                print(f"{project_id}. Extracting peri_events")
                count = 0
                for i in score_behaviors:
                    print(i)
                    behav=np.array(behav_scores[i]) 
                    if behav[start:end].any()==1:
                        store_peak=make_folder(i, peak_folder); store_behav=make_folder(i, behav_folder)
                        if i==lolo or i in point_events:
                            prestart=0
                            preend=fp_frames[-1]
                            placement=peri_event_analysis(behav, prestart, preend, i, l_base, behav_scores)
                            peri_baseline, peri_event, placement, placement_s=peri_event_splits(smooth_trace,placement)
                            event_means, baseline_means, peri_baseline_matrix, peri_event_matrix, placement=good_values(peri_baseline, peri_event, placement, i)
                            amplitude=[]; bigbad=[]; place_point=[]; probability=[]
                            save_me(peri_baseline_matrix, store_peak, project_id, peri_event_matrix, event_means, baseline_means, placement_s, i, amplitude, bigbad, place_point, probability)
                            save_me(peri_baseline_matrix, store_behav, project_id, peri_event_matrix, event_means, baseline_means, placement_s, i, amplitude, bigbad, place_point, probability)
                        if i!=lolo or i not in point_events:
                            placement=find_place(behav, start, end)
                            placement=peri_event_analysis(behav, start, end, i, l_base, behav_scores)
                            print(placement)
                            if len(placement)==0:
                                print(f'{i} is empty')
                            behav_analysis(placement, store_behav, project_id, behav_scores, score_behaviors, smooth_trace, fp_time, behav_time, fp_start, fp_end, i, animal_path, stage, "All Events",behav, count)
                            count+=1
                if use_groups=='yes':
                    make_groups(Groups, score_behaviors, behav_scores, start, end, l_base, peak_folder, behav_folder, heights, project_id, smooth_trace, fp_time, fp_start, fp_end, behav_time, animal_path, stage)
                
    print("Analysis complete. Subjects not analyzed: ", report)

def correct_photo_fps(_470, _410, behav_time, fp_time):

    #Find current frames per second for photometry data
    input_fps = int(len(_470)/behav_time)
    #Identify how many frames need to be added per second
    total_frames_togenerate = int(behav_time*fp_fps)

    new_index = np.linspace(0, len(_470) - 1, num=total_frames_togenerate, endpoint=True)
    #interpolate values
    interp_index = np.linspace(0, len(_470) - 1, total_frames_togenerate)
    new_470 = np.interp(interp_index, np.arange(len(_470)), _470)
    new_410 = np.interp(interp_index, np.arange(len(_410)), _410)

    return new_470, new_410, np.array([i/fp_fps for i in range(len(new_470))]), np.array([i for i in range(len(new_470))])
def make_groups(Groups, score_behaviors, behav_scores, start, end, l_base, peak_folder, behav_folder, heights, project_id, smooth_trace, fp_time, fp_start, fp_end, behav_time, animal_path,k):
    group_names=Groups.keys()
    group_lists=Groups.values()

    for list, name in zip(group_lists, group_names):    
        print(name)                  
        placement=[]
        for i in list:
            placement.extend(find_place(behav_scores[i],start, end))
        if len(placement)!=0:
            store_peak=make_folder(name, peak_folder); store_behav=make_folder(name, behav_folder)
            peak_analysis(placement, l_base, heights, store_peak, project_id, name, smooth_trace, behav_scores, score_behaviors, fp_time, behav_time, fp_start, fp_end, k, "Peaks", animal_path)
            placement=[]
            for i in list:
                placement.extend(peri_event_analysis(behav_scores[i],start, end, i, l_base, behav_scores))
            fix_place_2=placement[1:]; fix_place=placement[0:-1]
            for t in fix_place:
                for j in fix_place_2:
                    if int(j) in range(t+1, (t+(post_s*fp_fps))) or int(j-((pre_s-(floor(pre_s/2)))*fp_fps)) in range(t+1, (t+(post_s*fp_fps))):  #should prevent overlapping ROIs
                        if j in placement:
                            placement.remove(j)
            placement=[*set(placement)]

            behav = []
            behav_analysis(placement, store_behav, project_id, behav_scores, score_behaviors, smooth_trace, fp_time, behav_time, fp_start, fp_end, name, animal_path, k, "All Events", behav, 0)

def probable(l_base, maybe_mines, fp_start, fp_end, i, behav):

    start=int(floor(fp_start*behav_fps/fp_fps))
    end=(int(floor(fp_end*behav_fps/fp_fps)))
    
    placement=peri_event_analysis(behav, start, end, i, l_base, behav_scores)
    if len(maybe_mines)>1:
        z=len(maybe_mines)/len(placement) #peak probability probability that a peak will occur if that behavior happens   behaviors associated w peaks/all behaviors   
    else:
        z=np.nan
    if len(maybe_mines)<len(l_base):
        zz=len(maybe_mines)/len(l_base) #behav_probability  probability that the behavior will elicit a peak  behaviors associated w peaks/ALL peaks
    else:
        zz=1
    zzz=[z, zz]
    return zzz
def find_amp(smooth_trace, x,y):
    w=[smooth_trace[i] for i in x]
    z=[a-b for a, b in zip(y,w)]

    return z, w, y

def find_place(behav, start, end):
    placement=[]
    

    #IDENTIFIES THE BEGINNING OF EACH BEHAVIOR EVENT
    for i in range(len(behav)):
        check=int(i-1-(post_s*fp_fps))
        if start<i<end and check>0:
            a=i-1
            c=i-behav_fps
            b=i+floor(behav_fps/3)
            # c=i+(5*behav_fps)
            if behav[i]==1 and behav[a]==0 and behav[c]==0 and behav[i:b].any()==1:
                j=floor((i/behav_fps)*fp_fps)
                placement.append(j)   
            elif behav[i]==0:
                continue
    return placement

def behav_analysis(placement, store_project, project_id, behav_scores, score_behaviors, smooth_trace, fp_time, behav_time, fp_start, fp_end, i, animal_path, p, m, behav, count):
    
    peri_baseline, peri_event, placement, placement_s=peri_event_splits(smooth_trace,placement)
    event_means, baseline_means, peri_baseline_matrix, peri_event_matrix, placement=good_values(peri_baseline, peri_event, placement, i)
    amplitude=[]; bigbad=[]; place_point=[]; probability=[]; l_base=[]    
    pretty_picture(smooth_trace, fp_time, placement, fp_start, fp_end, i, m, l_base, animal_path,p, behav, count)
    save_me(peri_baseline_matrix, store_project, project_id, peri_event_matrix, event_means, baseline_means, placement_s, i, amplitude, bigbad, place_point, probability)
                    
def good_values(peri_baseline, peri_event, placement, i):

    event_means=[mean(peri_event[i]) for i in peri_event]
    baseline_means=[mean(peri_baseline[i][int(len(peri_baseline[i])/2):]) for i in peri_baseline]  #crops baseline for averaging
    
    diff=[(100/a*b-100) for a,b in zip(baseline_means, event_means)]
    if stop_at>0:     
        try:
            workme=pd.DataFrame({'event':event_means, 'baseline':baseline_means, 'placement':placement,'base2': list(peri_baseline), 'event2': list(peri_event), 'diff': diff})
            workme2=workme.sort_values(by='diff', ascending=[False])
            if i in Groups:
                df=workme2[0:int(stop_at*2)]
            else:
                df=workme2[0:int(stop_at)]
            
            event_means=df['event']
            baseline_means=df['baseline']
            placement=df['placement']
            
            peri_baseline_matrix=np.array([peri_baseline[i] for i in df['base2']], dtype=float)  #keeps equal length values for AUC
            peri_event_matrix=np.array([peri_event[i] for i in df['event2']], dtype= float)  #keeps equal length values for AUC
        except ValueError:
            event_means=[]; baseline_means=[]; peri_baseline_matrix=[]; peri_event_matrix=[]; placement=[]
    else:
        peri_baseline_matrix=np.array([peri_baseline[i] for i in peri_baseline], dtype=float)  #keeps equal length values for AUC
        peri_event_matrix=np.array([peri_event[i] for i in peri_event], dtype= float)  #keeps equal length values for AUC
        

    return event_means, baseline_means, peri_baseline_matrix, peri_event_matrix, placement;                         
                        
def zones_in_use(sniff_group, behav_raw, score_behaviors, zone):

    behav_scores={}
    work=[]
    #MAKE AS BOOLEAN ARRAY
    for t in score_behaviors:
        work.append(t)
        j=np.array(behav_raw[t], dtype=bool)  ##crops behavior by cutoff time
        behav_scores.update({f"{t}": j})
    for s in sniff_group:
        yes=[0] * fp_fps
        no=[0] * fp_fps
        for p,q in zip(behav_raw[s], behav_raw[zone]):
            if p==1 and q==1:
                yes.append(1)
            else:
                yes.append(0)
            if p==1 and q==0:
                no.append(1)
            if p==1 and q==1:
                no.append(0)
            else:
                no.append(0)
        j=np.array(yes, dtype=int)
        k=np.array(no, dtype=int)
        l=str(f"{s}_inZone")
        print(l)
        m=str(f"{s}_OutZone")
        work.append(l); work.append(m)
        behav_scores.update({f"{l}": j})
        behav_scores.update({f"{m}": k})
        
        print(score_behaviors)

    return behav_scores, work;

def save_me(peri_baseline_matrix, store_project, project_id, peri_event_matrix, event_means, baseline_means, placement_s, i, amplitude, bigbad, place_point, probability):
    if i!=lolo or i not in point_events:
        if count(baseline_means)!=0 and count(event_means)!=np.nan and count(baseline_means)!=np.nan:
            means_to_compare=[mean(baseline_means), mean(event_means)]
            if count(means_to_compare)!=0 or means_to_compare!=np.nan: 
                np.savetxt(f'{store_project}/means_to_compare_{project_id}.csv', means_to_compare, delimiter=',', fmt='%s')
        if len(amplitude)!=0:
            np.savetxt(f'{store_project}/amplitude_{project_id}.csv', amplitude, delimiter=',', fmt='%s')
        if len(bigbad)!=0:
            np.savetxt(f'{store_project}/peak_prom_{project_id}.csv', bigbad, delimiter=',', fmt='%s')
        if len(place_point)!=0:
            np.savetxt(f'{store_project}/place_value_{project_id}.csv', place_point, delimiter=',', fmt='%s')
        if len(probability)!=0:    
            np.savetxt(f'{store_project}/probability_{project_id}.csv', probability)

    if count(peri_baseline_matrix)!=0:
        np.savetxt(f'{store_project}/peri_baseline_matrix_{project_id}.csv', peri_baseline_matrix, delimiter=',', fmt= '%s')
    if count(peri_event_matrix)!=0:
        np.savetxt(f'{store_project}/peri_event_matrix_{project_id}.csv', peri_event_matrix, delimiter=',', fmt='%s')  #These two will have equal length values for AUC and z-score calculations
    if count(event_means)!=0:  #starting here, baseline is cropped to 'pre_s' length
        np.savetxt(f'{store_project}/event_means_{project_id}.csv', event_means, delimiter=',', fmt='%s')
    if count(baseline_means)!=0:
        np.savetxt(f'{store_project}/baseline_means_{project_id}.csv', baseline_means, delimiter=',', fmt='%s')
    if count(placement_s)!=0: 
        np.savetxt(f'{store_project}/placement_s_{project_id}.csv', placement_s, delimiter=',', fmt='%s')

def plot_me(behav_scores, score_behaviors, smooth_trace, fp_time, behav_time, fp_start, fp_end, l_base, maybe_mines, i, animal_path, p, m):

    data = np.array([behav_scores[w] for w in score_behaviors if w!=lolo])
    x_axis = [m for m in range(len(behav_scores[lolo]), 10*fp_fps)]
    x_labels = [m/fp_fps for m in x_axis]

    # create some labels
    labels = [ v for v in score_behaviors]
    # labels = ["Tail suspension","Rear","Hand chase", " "]
    color_list=['lawngreen', 'dodgerblue', 'orangered', 'navy', 'magenta', 'green', 'gold','deeppink', 'aquamarine','darkviolet', 'olivedrab', 'deepskyblue','palevioletred', 'orange', 'mediumturquoise', 'red','sandybrown', 'orchid', 'springgreen', 'darkorchid', 'royalblue']

    # # create a color map with random colors
    colmap = matplotlib.colors.ListedColormap(color_list)
    
    colmap.colors[0] = [1,1,1]

    # # create some colorful data:
    data_color = (1 + np.arange(data.shape[0]))[:, None] * data

    y_=len(labels)

    centers=[behav_time.min(), behav_time.max(), y_+1,0]
    dx, = np.diff(centers[:2])/(data_color.shape[1]-1)
    dy, = -np.diff(centers[2:])/(data_color.shape[0]-1)

    extent = [centers[0], centers[1], centers[2]-dy, centers[3]-dy/2]

    #ETHOGRAM ALIGNED WITH NORMALIZED TRACE
    etho_fig = plt.figure()
    gs = etho_fig.add_gridspec(2, hspace=0)
    ax = gs.subplots(sharex=True)

    # x_label = list(np.arange(0,len(behav_scores[lolo]),(10*behav_time)))

    # x_labels = [x/fp_fps for x in x_label]

    ax[1].imshow(data_color, aspect='auto',cmap=colmap, interpolation='nearest', extent = extent)
    ax[1].set_yticks(np.arange(len(labels)))
    ax[1].set_yticklabels(labels)  ##fix y labels and tighten the graph
    # ax[1].set_xticks(range(0, len(behav_scores[lolo]), 10*fp_fps))
    ax[1].set_xticklabels(x_labels)
    # plt.xticks(x_axis, x_labels)
    etho_fig.align_ylabels()


    ax[0].plot(fp_time,smooth_trace)
    # ax[0].label_outer()
    ax[0].set_ylim([smooth_trace.min()-(0.30*smooth_trace.min()), smooth_trace.max()+(0.30*smooth_trace.max())])
    
    ax[0].set_title(f'{i}')

    x_positions = np.arange(0,fp_time.max(),(100*fp_fps))
    # etho_fig.xticks(x_labels)  

    plt.tight_layout()


        # print(placement_s)
    for k in maybe_mines:
        a=int(k)
        b=int(k+post_s*20)
        c=int(k-pre_s*20)
        ax[0].axvspan(fp_time[a], fp_time[b], color='r',  alpha=0.4, lw=0)
        ax[0].axvspan(fp_time[a], fp_time[c], color='b', alpha=0.4, lw=0)
    if m=="Peaks":
        ax[0].scatter(fp_time[l_base], smooth_trace[l_base], color='g', s=30)
    ax[0].axvspan(fp_time[fp_start], fp_time[fp_end], color='y', alpha=0.2, lw=0)
    # ax[0].scatter(fp_time[r_base], smooth_trace[r_base], color='b', s=30)
    # ax[0].plot(fp_time[prom], smooth_trace[prom], color='r', s=30)
    # if p=="during": 
    # plt.show()
    roi_fig = make_folder("ROI_figures", animal_path)
    roi_fig2=make_folder(m, roi_fig)

    etho_fig.savefig(f'{roi_fig2}\{i}_{p}_ROI.tif')

    plt.close('all')

def pretty_picture(smooth_trace, fp_time, maybe_mines, fp_start, fp_end, i, m, l_base, animal_path,p, behav, count):
    for_plot_start = []
    for_plot_end = []

    for j in range(len(behav)):
        if behav[j-1]==False and behav[j]==True:
            for_plot_start.append(round(j/behav_fps*fp_fps))
        if behav[j]==True and behav[j+1]==False:
            for_plot_end.append(round(j/behav_fps*fp_fps))
    
    
    # color_list=['green', 'gold','red', 'lawngreen', 'dodgerblue',  'green', 'deeppink', 'aquamarine', 'deepskyblue','palevioletred', 'orange', 'mediumturquoise', 'red','sandybrown', 'orchid', 'springgreen', 'darkorchid', 'royalblue']
    color_list = ['dodgerblue', 'lawngreen', 'gold', 'magenta']
    # # create a color map with random colors
    colmap = matplotlib.colors.ListedColormap(color_list)


    #ETHOGRAM ALIGNED WITH NORMALIZED TRACE
    plt.figure(figsize = (10,5))

    # plt.align_ylabels()


    plt.plot(fp_time,smooth_trace, color = 'k')
    # etho_fig.label_outer()
    plt.ylim([smooth_trace.min()-(0.30*smooth_trace.max()), smooth_trace.max()+(0.10*smooth_trace.max())])
    plt.ylabel(r'zF', fontsize = 20)
    plt.xlabel('Time (s)', fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    # plt.set_title(f'{i}')

    x_positions = np.arange(0,fp_time.max(),(100*fp_fps))
    # etho_fig.xticks(x_labels)  
    # print(new_behav_dict)
    plt.tight_layout()

    plt.axvspan(fp_time[fp_start], fp_time[fp_end], color='palegoldenrod', alpha=0.2, lw=0)
    for k in maybe_mines:
        a=int(k)
        b=int(k+post_s*fp_fps)
        c=int(k-pre_s*fp_fps)
        plt.axvspan(fp_time[a], fp_time[b], color='coral',  alpha=0.5, ymax = 0.6, label = i)
        plt.axvspan(fp_time[a], fp_time[c], color='aquamarine',  alpha=0.5, ymax = 0.6, label = i)
    if m=="Peaks":
        plt.scatter(fp_time[l_base], smooth_trace[l_base], color='g', s=30)

    if len(behav)!= 0:
        for f, h in zip(for_plot_start, for_plot_end):
            try:
                plt.axvspan(fp_time[f], fp_time[h], color='midnightblue', ymax = 0.1)
            except IndexError:
                continue

    plt.title(f"{i}")
    roi_fig = make_folder("ROI_figures", animal_path)
    roi_fig2=make_folder(m, roi_fig)

    plt.savefig(f'{roi_fig2}\{i}_{p}_ROI.tif')
    # plt.show()
    plt.close('all')


    # plt.show()
    # plt.close('all')
def peri_event_splits(smooth_trace, placement):
  
    peri_baseline={}
    peri_event={}

    #STORES THE BASELINE AND PERI-EVENT
    counter=1
    for i in placement:
        m=int(i-1-(post_s*fp_fps))
        n=int(i-1)
        q = int(i+(post_s*fp_fps))
        if m>0:
            baseline=[]
            event=[]
            base_start=[]
        #RETRIEVES VALUES FROM THE SMOOTH NORMALIZED TRACE
            for j, k in zip(smooth_trace[m:n], smooth_trace[i:q]): #Overshoot for AUC
                if m >= 0 and q+5 <= len(smooth_trace):
                    baseline.append(j)
                    event.append(k)
            peri_baseline.update({f'{counter}':baseline})
            peri_event.update({f'{counter}': event})
            counter+=1
    placement_s=[floor((i)/fp_fps) for i in placement]

    return peri_baseline, peri_event, placement, placement_s;

def peri_event_analysis(behav, start, end, i, l_base, behav_scores):
    placement=[]

    #IDENTIFIES THE BEGINNING OF EACH BEHAVIOR EVENT

    for o in range(len(behav)):
        if start<o<end:
            a=o-1
            b=int(o+floor((duration_s)*behav_fps))
            c=o+1

            if i not in point_events:
                if behav[o]==1 and behav[a]==0 and behav[o:b].all()==1:
                    j=floor((o/behav_fps)*fp_fps)
                    placement.append(j-offset*fp_fps)   
                if behav[o]==1 and behav[c]==1:
                    # print('yes')
                    j=floor((o/behav_fps)*fp_fps)
                    if j in l_base: # and j not in check_first:
                        placement.append(j-offset*fp_fps)  
            else:
                if behav[o]==1 and behav[a]==0:
                    j=floor((o/behav_fps)*fp_fps)
                    placement.append(j-offset*fp_fps) 
                elif behav[o]==0:
                    continue

    fix_place_2=placement[1:]; fix_place=placement[0:-2]

    for d,j in zip(fix_place, fix_place_2):
        if int(j) in range(d, (d+(pre_s*fp_fps))) or int(j-((pre_s-(floor(pre_s/2)))*fp_fps)) in range(d, (d+(floor(post_s/1.5)*fp_fps))):  #should prevent overlapping ROIs
            if j in placement:
                placement.remove(j)

    return placement

def peak_analysis(placement, l_base, heights, store_project, project_id, i, smooth_trace, behav_scores, score_behaviors, fp_time, behav_time, fp_start, fp_end, p, m, animal_path):

    amplitude=[]; bigbad=[]; place_point=[]; probability=[]
    maybe_mines=[] #new placements
    bigbad=[] #new heights

    for a, b in zip(l_base, heights):
        for p in placement:
            check=int(p-1-(post_s*fp_fps))
            if p in range((a-(floor(pre_s)*fp_fps)),a+(floor(pre_s)*fp_fps)) and p not in maybe_mines and check>0:# and p not in check_first:
                maybe_mines.append(p)

                bigbad.append(b)
                break
            else:
                continue
    
    fix_place_2=maybe_mines[1:]; fix_place=maybe_mines[0:-1]

    for h in fix_place:
        for j in fix_place_2:
            pop=1
            if int(j) in range(h+2, (h+(post_s*fp_fps))) or floor(j-((pre_s-(floor(pre_s/2)))*fp_fps)) in range(h+2, (h+(post_s*fp_fps))):  #should prevent overlapping ROIs
                if j in maybe_mines:
                    maybe_mines.remove(j)
                    bigbad.remove(bigbad[pop])
            pop+=1

    if len(placement)>0:
        # probability=probable(l_base, maybe_mines, fp_start, fp_end, i, behav)
        peri_baseline, peri_event, placement, placement_s=peri_event_splits(smooth_trace,maybe_mines)
        event_means, baseline_means, peri_baseline_matrix, peri_event_matrix, placement=good_values(peri_baseline, peri_event, placement, i)
        # plot_me(behav_scores, score_behaviors, smooth_trace, fp_time, behav_time, fp_start, fp_end, l_base, placement, i, animal_path, p, m)
        amplitude, place_point, bigbad=find_amp(smooth_trace, placement, bigbad)
        save_me(peri_baseline_matrix, store_project, project_id, peri_event_matrix, event_means, baseline_means, placement_s, i, amplitude, bigbad, place_point, probability)

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

def cropendtime(_470, _410, behav_raw, behavior, crop_end):
    _470_=np.array([_470[i] for i in range(len(_470)) if i<(len(_470)-(crop_end*fp_fps))])
    _410_=np.array([_410[i] for i in range(len(_410)) if i<(len(_410)-(crop_end*fp_fps))])
    frame = np.array([i for i in range(len(_470_))])
    time = np.array([i/fp_fps for i in frame])

    behav_raw_=behav_raw[:-(crop_end*behav_fps)]
    behav_frames = [i for i in range(len(behav_raw_[behavior[0]]))]
    # TIME
    behav_time = np.array([i/behav_fps for i in behav_frames])

    return _470_,_410_, behav_raw_, behav_time, time, frame, behav_frames;

def crop_trace(_470, _410, behav_raw, behavior, start, end):
    fp_start=int(start/behav_fps*fp_fps); fp_end=int(end/behav_fps*fp_fps)
    mini=start-time_from_start*behav_fps
    end_length = len(_470)
    if mini>0:
        behav_raw_=behav_raw[start-time_from_start*behav_fps: end+time_from_end*behav_fps]
        new_470=_470[fp_start-time_from_start*fp_fps:fp_end+time_from_end*fp_fps]
        new_410=_410[fp_start-time_from_start*fp_fps:fp_end+time_from_end*fp_fps]
    elif mini<0:
        # mini2=0-mini
        # mini3=0-(fp_start-time_from_start*fp_fps)
        behav_raw_=behav_raw[:end+time_from_end*behav_fps]
        new_470=_470[:fp_end+time_from_end*fp_fps]
        new_410=_410[:fp_end+time_from_end*fp_fps]
    elif end_length<(fp_end+(time_from_end*fp_fps)):
        behav_raw_=behav_raw[start-time_from_start*behav_fps:]
        new_470=_470[fp_start-time_from_start*fp_fps:]
        new_410=_410[fp_start-time_from_start*fp_fps:]

    # _470_=np.array([i for i in new_470])
    # _410_=np.array([i for i in new_410])
    frame = np.array([i for i in range(len(new_470))])
    time = np.array([i/fp_fps for i in frame])

    behav_frames = [i for i in range(len(behav_raw_[behavior[0]]))]
    # TIME
    behav_time = np.array([i/behav_fps for i in behav_frames])

    return new_470,new_410, behav_raw_, behav_time, time, frame, behav_frames;

def takeoff(behav_time, fp_time, behavior, behav_fps, behav_raw):
    
    # fp_addon_offset=behav_duration-fp_time[-1]
    # #rationale: the behavior starts before photometry; but the behavior video decides when the photometry should start
    cutit = 0
    if behav_time.max()+2>fp_time.max():
        offset=(behav_time.max()-fp_time.max())
        cutit=floor(offset*behav_fps) ##Make sure it is rounding down

        j=behav_raw[behavior[0]]
        behav_frames=np.array([i for i in range(len(j[cutit:]))])
        behav_times=np.array([i/behav_fps for i in behav_frames])
        behav_raw=behav_raw[cutit:]
        behav_time = behav_times

    return behav_time, behav_raw, cutit;

def behav_split(behav_raw, behav_fps):
    # GET BEHAVIOR HEADERS
    behaviors = [i for i in behav_raw]
    # only registers columns 8 to the second to last as behaviors; ignores 'Result 1'
    behavior = behaviors[7:-1]

    # FRAME NUMBERS
    behav_frames = [i for i in range(len(behav_raw[behavior[0]]))]

    # TIME
    behav_time = np.array([i/behav_fps for i in behav_frames])

    return behavior, behav_frames, behav_time

def photo_info(_470, _410, fp_frames, project_id):

    # CURVE FITTING
    fit1, fit2, fit3 = fit_that_curve(fp_frames, _410, _470)

    # CORRECT FOR POSSIBLE LED STATE SWTICH  ##doesnt work if 410 is strong
    _410_ = (_410[1000:])*(1/_410[1000].max())
    _470_ = _470[1000:]*(1/_410[1000].max())
    diff1 = _410_.max()-_410_.min()
    diff2 = _470_.max()-_470_.min()
    counts = 0

    for i in _410_:
        if i == _410_.max():
            counts = +1
    # print(count)

    if diff1 > diff2 or counts > 1:
        _410_ = _470
        _470 = _410
        _410 = _410_

        fit1, fit2, fit3 = fit_that_curve(fp_frames, _410, _470)

    # NORMALIZE AND SMOOTH TRACE
    normalized = (_470-fit3)/fit3
    smooth_normal = smooth(normalized, 59)  # must be an odd number
    smooth_trace_ = smooth_normal-smooth_normal.min()
    scale = 1/smooth_trace_.max()
    smooth_trace = smooth_trace_*scale

    plot_photo(fp_frames, _470, _410, project_id, fit1, fit3, smooth_trace)
# 
    while True:
        correct=input("Is 410 and 470 correct? (Y/N) Or press (S) to skip: ")
        # correct = 'y'
        plt.close()
        if correct.lower()=='y' or correct.lower()=='':
            break
        if correct.lower()=='n':
            _410_=_470; _470=_410; _410=_410_
            print("switching 410 and 470")
            break
        if correct.lower()=='s':
            _ = []; __ = []; ___ = []
            return correct, _, __, ___;
    fit1, fit2, fit3= fit_that_curve(fp_frames, _410, _470)
    # NORMALIZE AND SMOOTH TRACE
    while True:
        good_control = input("Remove isobesctic control? (y/n): ")  
        if good_control.lower() == 'y':
            fp_times=[i/fp_fps for i in fp_frames]
            normalized=(_470-fit3)/fit3
            smooth_normal = smooth(_470, 59)  # must be an odd number
            smooth_trace_ = smooth_normal-smooth_normal.min()
            scale=1/smooth_trace_.max()
            smooth_trace=smooth_trace_*scale
            break
        if good_control.lower() == 'n' or good_control.lower() =='':
            fp_times=[i/fp_fps for i in fp_frames]
            normalized=(_470-fit3)/fit3
            smooth_normal = smooth(normalized, 59)  # must be an odd number
            smooth_trace_ = smooth_normal-smooth_normal.min()
            scale=1/smooth_trace_.max()
            smooth_trace=smooth_trace_*scale
            break

    # plot_photo(fp_frames, _470, _410, project_id, fit1, fit3, smooth_trace)

    # LIST OF PEAK POSITIONS
    x, _ = list(signal.find_peaks(smooth_trace, width=50, height=0.21))
    prom, l_base, r_base = list(signal.peak_prominences(smooth_trace, x, wlen=300))
    heights = smooth_trace[x]

    return smooth_trace, fp_frames, heights, l_base;

def open_files(file_to_load, animal_path, animal, report):
    if len(file_to_load)>=2:
        for i in file_to_load:
            try:
                if i.endswith('.csv') or i.endswith('.CSV'):
                    file = rf'{animal_path}\{i}'
                    for i in file:
                        if i == '\\':
                            file.replace(i, '\\\\')
                    fp_raw = pd.read_csv(file)
                    # files.append(file)
                    continue
                if 'Raw data' in i and i.endswith('.xlsx') and '$' not in i:
                    file = rf'{animal_path}\{i}'
                    for i in file:
                        if i == '\\':
                            file.replace(i, '\\\\')
                    behav_raw = pd.read_excel(file, header=[behavior_row-1], skiprows=[behavior_row])
                    # behav_raw=pd.read_excel(file, header=[34], skiprows=[35])
                    behavior, behav_frames, behav_time = behav_split(
                        behav_raw, behav_fps)

                    if lolo not in behavior:
                        new_row=int(input("Behaviors are not in designated row. Try again: "))
                        behav_raw = pd.read_excel(file, header=[new_row-1], skiprows=[new_row])
                        # behavior_row=new_row
            except ValueError:
                print(f"Unable to open file: {file}")
                report.append(animal)
                continue
    return fp_raw, behav_raw, report;


def open_gui(example_behav_file, behavior_row):
    # GUI FOR BEHAVIOR SELECTION
    
    while True:
        question=input("Load or Select behaviors? (Respond L for Load, S for Select): ")
        # question = 'l'
        if question.lower() == 'l':
            with open('select_behaviors.csv', 'r') as read_obj:
                score = list(csv.reader(read_obj))
            score_behaviors = [i for i in score]
            score_behaviors = score_behaviors[0]
            print(score_behaviors)
            # score_behaviors
            break
        elif question.lower() == 's':
            find = pd.read_excel(example_behav_file, header=[behavior_row-1], skiprows=[behavior_row])
            behavior_exp, behav_frames_exp, behav_time_exp = behav_split(
                find, behav_fps)
            # print(behavior_exp)
            if lolo not in behavior_exp:
                find = pd.read_excel(example_behav_file, header=[behavior_row-2], skiprows=[behavior_row-1])
                if lolo not in behavior_exp:
                    find = pd.read_excel(example_behav_file, header=[behavior_row], skiprows=[behavior_row+1])
                    if lolo not in behavior_exp:
                        new_row=int(input("Behaviors are not in designated row. Try again: "))
                        find = pd.read_excel(example_behav_file, header=[new_row-1], skiprows=[new_row])
                behavior_exp, behav_frames_exp, behav_time_exp = behav_split(
                    find, behav_fps)
                behavior_row=new_row
            score_behaviors = pick_behaviors(behavior_exp)  # GUI
            # score_behaviors
            break
        else:
            continue
    return score_behaviors

def fit_that_curve(fp_frames, _410, _470):
    # CALCULATIONS FOR BIEXPONENTIAL FIT
    popt, pcov = curve_fit(exp2, fp_frames, _410, maxfev=500000, p0=[
                           0.01, 1e-7, 1e-5, 1e-9])
    fit1 = exp2(fp_frames, *popt)

    # LINEARLY SCALE FIT TO 470 DATA USING ROBUST FIT
    # https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
    A = np.vstack([fit1, np.ones(len(fit1))]).T
    # https://www.mathworks.com/help/stats/robustfit.html
    slope = np.linalg.lstsq(A, _470, rcond=None)[0]
    fit2 = fit1*slope[0]+slope[1]

    # FIT 410 OVER 470
    fit3_ = stats.linregress(_410, _470)
    fit3 = fit3_.intercept+fit3_.slope*_410

    return fit1, fit2, fit3

def exp2(x, a, b, c, d):
    return a*exp(b*x) + c*exp(d*x)

# GET PHOTOMETRY VALUES FROM EXCEL SHEET
def fp_split(fp_raw, z, fp_fps):

    roi = []
    # Extract column headers from CSV to give user options
    for i in fp_raw:
        if 'Region' in i:
            roi.append(i)
        else:
            continue
    # Indexing column headers, not values
    roi_ = roi[z]

    led = fp_raw['LedState']
    trace = fp_raw[roi_]

    # print(trace)
    _470_ = []
    _410_ = []
    c = 0
    if full_trace=='yes':
        for i, j in zip(led, trace):
            if i == 6 and c > (crop_front*fp_fps):
                _470_.append(j)
            elif i != 6 and c > (crop_front*fp_fps):
                _410_.append(j)
            c += 1
    else:
        for i, j in zip(led, trace):
            if i == 6:
                _470_.append(j)
            elif i != 6:
                _410_.append(j)
            c += 1

    if len(_470_) != len(_410_):
        if len(_470_) > len(_410_):
            _470_.remove(_470_[-1])
        elif len(_470_) < len(_410_):
            _410_.remove(_410_[-1])

    _470 = np.array(_470_)
    _410 = np.array(_410_)

    frame = np.array([i for i in range(len(_470))])

    time = np.array([i/fp_fps for i in frame])

    return time, _470, _410

def plot_photo(fp_frames, _470, _410, project_id, fit1, fit3, smooth_trace):

    # #plotS
    # RAW TRACE
    fig, ax = plt.subplots(4)
    ax[0].plot(fp_frames, _470)
    ax[0].set_title('Raw 47\0')
    ax[0].set(xlabel='Frame #', ylabel=r'$\Delta$F/F')

    # BIEXPONENTIAL FIT
    ax[1].set_ylim([_410.min(), _410.max()])
    ax[1].plot(fp_frames, _410)
    ax[1].plot(fp_frames, fit1, 'r')
    # rationale: https://onlinelibrary.wiley.com/doi/10.1002/mma.7396
    ax[1].set_title('410 with biexponential fit')
    ax[1].set(xlabel='Frame #', ylabel=r'$\Delta$F/F')
    ax[0].set_title(f'{project_id}')
    # ax[3].set(xlabel='Frame #', ylabel=r'F mean pixels')

    # 410 FIT TO 470
    ax[2].plot(fp_frames, _470)
    ax[2].plot(fp_frames, fit3, 'r')
    ax[2].set_title('410 fit over 470')
    ax[2].set(xlabel='Frame #', ylabel=r'F mean pixels')

    ax[3].plot(fp_frames, smooth_trace)
    plt.show(block=False)

def pick_behaviors(behavior):  # from Tkinter import * #https://www.geeksforgeeks.org/how-to-create-a-pop-up-message-when-a-button-is-pressed-in-python-tkinter/
    root = tkinter.Tk()
    root.title("Select which behaviors to calculate peri-events")
    root.geometry('500x300')

    # https://stackoverflow.com/questions/50485891/how-to-get-and-save-checkboxes-names-into-a-list-with-tkinter-python-3-6-5
    class CheckBox(tkinter.Checkbutton):
        boxes = []  # Storage for all buttons

        def __init__(self, master=None, **options):
            # Subclass checkbutton to keep other mbehavds
            tkinter.Checkbutton.__init__(self, master, options)
            self.boxes.append(self)
            self.var = tkinter.BooleanVar()  # var used to store checkbox state (on/off)
            self.text = self.cget('text')  # store the text for later
            # set the checkbox to use our var
            self.configure(variable=self.var)
    for i in behavior:
        c = CheckBox(root, text=i).pack()

    save_button = tkinter.Button(root, text="SAVE", command=root.destroy)
    save_button.pack()

    # EXECUTES GUI
    root.mainloop()
    score_behaviors = []
    # SAVES SELECTIONS
    for box in CheckBox.boxes:
        if box.var.get():  # Checks if the button is ticked
            score_behaviors.append(box.text)

    with open('select_behaviors.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(score_behaviors)
        
    return score_behaviors

def count(j):
    k = 0
    for i in range(len(j)):
        k += 1
    return k

# SMOOTH AND NORMALIZE PHOTOMETRY TRACE
def smooth(a, WSZ):  # https://stackoverflow.com/questions/40443020/matlabs-smooth-implementation-n-point-moving-average-in-numpy-python
    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid')/WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((start, out0, stop))

if __name__ == "__main__":  # avoids running main if this file is imported
    main()