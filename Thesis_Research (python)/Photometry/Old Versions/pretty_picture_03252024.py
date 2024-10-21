# AUTHOR: VICTORIA SEDWICK
# ADAPTED FROM NEUROfpMETRICS AND ILARIA CARTA (MATLAB)

# https://mat#plotlib.org/stable/gallery/sub#plots_axes_and_figures/sub#plots_demo.html
# https://stackoverflow.com/questions/11190735/python-mat#plotlib-superimpose-scatter-#plots

# <Where are your files located?>


fp_folder = r"D:\Tests\82"
# project_home=input("Project folder?: ")
# NOTE <Load any ethovision file that has behavior names listed>
# example_behav_file = r"E:\Photometry-Fall2022\Pup Block\Virgin\Master file (Males)\Infanticidal Animals_trials\Pups\10\Raw data-MeP-CRFR2photo_pupblock1-Trial    31 (1).xlsx"
# <Frames per second for photo and ethovision data>
behav_fps = 30
fp_fps = 20
# <How many seconds do you want to analyze for peri-event?>
pre_s = 5
post_s = 10

crop_end=10; crop_front=10
duration_s=0.5; behavior_row=33
#YES OR NOl

full_trace='no'  #if no, the program will crop around the start and end times
use_zones='no'
use_groups='no'
controls='no'

#only if full_trace=no
time_from_start=60; time_from_end=30
    
stop_at=0

# Which roi/fiber signal do you want to #plot (0-2)?
z = 0
# name of start behavior
lolo = str('start')
first = str('Approach')
point_events=['Approach', 'start', 'bite']

score_behaviors = ['Attack','Aggressive groom','Groom','Sniff', 'Retrieval', 'start', 'Carry']

colored_behaviors = {'Retrieval': 'magenta', 
                     'Sniff': 'dodgerblue', 
                     'Groom': 'green', 
                     'Attack': 'red', 
                     'Aggressive groom': 'orangered', 
                     'Carry': 'aquamarine'}

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
import matplotlib.patches as mpatches


def main():
    global colored_behaviors, score_behaviors, fp_folder, behav_fps, behavior_row, fp_fps, example_behav_file, z, pre_s, post_s, lolo, controls, Groups, use_zones, use_groups, zone, crop_front, duration_s, min_duration, point_events, first, contact_behaviors, stop_at, complex_behaviors;


    # LOAD EXAMPLE BEHAVIOR FILE FOR LIST OF AVAILABLE BEHAVIORS
    report=[]
    rootdir = os.listdir(fp_folder)
    fp_file = ''; behav_file = '';
    for i in rootdir:
        ii = i.lower()
        if ii.endswith('.csv'):
            fp_file+=os.path.join(fp_folder, i)
        if ii.endswith('.xlsx') and '$' not in i:
            behav_file+=os.path.join(fp_folder, i)
    try:
        fp_raw, behav_raw = open_files(fp_file, behav_file)
    except UnboundLocalError:
        exit()

    fp_time, _470, _410 = fp_split(fp_raw, z, fp_fps)  #cuts the first 30 seconds
    behavior, behav_frames, behav_time = behav_split(behav_raw, behav_fps)

    print(f"Time difference between video and trace before cropping: {behav_time.max()}, {fp_time.max()}")

    #ALIGN BEHAVIOR AND PHOTOMETRY
    behav_time, behav_raw, cut=takeoff(behav_time, fp_time, behavior, behav_fps, behav_raw)  #wanna crop end time (100?)
    #START TIMES & CONTROLS OR NO CONTROLS
    start_times=np.array(behav_raw[lolo])
    place=[int(i) for i in range(len(start_times)) if start_times[i]==1]
    
    if full_trace.lower()=='no':
        print('yes')
        _470,_410, behav_raw, behav_time, fp_time, fp_frames, behav_frames=crop_trace(_470, _410, behav_raw, behavior, place[0], place[-1])
    else:           
        _470,_410, behav_raw, behav_time, fp_time, fp_frames, behav_frames=cropendtime(_470, _410, behav_raw, behavior, crop_end) #cuts the last minute

    #START TIMES & CONTROLS OR NO CONTROLS
    start_times=np.array(behav_raw[lolo])
    place=[int(i) for i in range(len(start_times)) if start_times[i]==1]

    start_placements = [1]
    end_placements = [place[0]]
    stages = ['before', 'during', 'controls']
    
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

    print(f"Time difference between video and trace after cropping: {behav_time.max()}, {fp_time.max()}")

    print('Normalizing photometry trace')
    smooth_trace, scaled_trace, zscore_trace, fp_frames, heights, l_base = photo_info(_470, _410, fp_frames)       
   
    trace_diary = {'Scaled_trace': scaled_trace, 'Zscore_trace': zscore_trace, 'Normalized trace': smooth_trace}
    fp_behav_scores = {}
    for i in score_behaviors:
        placements = [int(p/behav_fps*fp_fps) for p,_ in zip(range(len(behav_raw[i])), behav_raw[i]) if _ == 1]
        fp_array = []
        for frame in range(len(fp_frames)):
            if frame in placements:
                fp_array.append(1)
            else:
                fp_array.append(0)
        fp_behav_scores.update({f'{i}': np.array(fp_array, dtype = bool)})

    use_place = [use for use in zip(start_placements,end_placements,stages) if 'during' in use]   
    start, end, stage = use_place[0]
    fp_start = int(start*fp_fps/behav_fps)
    fp_end = int(end*fp_fps/behav_fps)
    label_directory = fp_folder.split('\\')
    assay_folder = '\\'.join(label_directory[:-3])
    print(label_directory[-3:])
    assay, subject, id_ = label_directory[-3:]
    Representative_image = make_folder("Representative Trace", "D:/")

    # print(behav_scores)
    for trace in trace_diary:
        plot_ethogram_trace(fp_behav_scores, trace_diary[trace], fp_time, fp_start, fp_end, subject, id_, Representative_image, trace, assay)


    
    print("Analysis complete. Subjects not analyzed: ", report)



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
                j=int((i/behav_fps)*fp_fps)
                placement.append(j)   
            elif behav[i]==0:
                continue
    return placement

def behav_analysis(placement, behav_scores, score_behaviors, smooth_trace, fp_time, behav_time, fp_start, fp_end, i):

    peri_baseline, peri_event, placement, placement_s=peri_event_splits(smooth_trace,placement)
    new_array = []
    for i in range(len(fp_time)):
        if i in placement:
            new_array.append(1)
        else:
            new_array.append(0)
    return new_array

def plot_ethogram_trace(new_behav_dict, trace, fp_time, fp_start, fp_end, subject, id_, Representative_image,  trace_name, assay):
    # Set up the figure and axis for the plot with gridspec_kw to adjust relative heights
    fig, (ax_trace, ax_ethogram) = plt.subplots(2, 1, figsize=(8, 4), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    plt.subplots_adjust(hspace = 0.1)
    ax_trace.tick_params(axis='y', labelsize=12)
    
    # Plot the photometry trace with a skinnier line
    ax_trace.plot(fp_time, trace, 'k', linewidth=0.5)  # 'k' for black line
    ax_trace.axvline(fp_time[fp_start], color = 'k', linestyle = '--')
    # ax_trace.set_title(f'{subject} - {assay}')
    
    if trace_name == 'Scaled_trace' or trace_name =='Normalized trace':
        ax_trace.set_ylabel(r'$\Delta$F', size = 15)
    elif trace_name == 'Zscore_trace':
        ax_trace.set_ylabel('zF', size = 15)
    
    # Set the y-limits to be tighter around the data
    trace_min, trace_max = min(trace), max(trace)
    ax_trace.set_ylim([trace_min - 0.1 * abs(trace_min), trace_max + 0.1 * abs(trace_max)])
    
    # Plot the ethogram below the photometry trace
    ethogram_base = -0.1  # Start plotting ethogram events at this y-value
    linelengths = 0.1  # Make the event lines shorter
    count = 0
    legend_patches = []
    for i, behavior in enumerate(sorted(colored_behaviors.keys())):
        behavior_data = new_behav_dict.get(behavior, [])
        if any(behavior_data):
            # Extract the times at which the behavior occurs
            behavior_times = fp_time[behavior_data]
            ax_ethogram.eventplot(behavior_times, lineoffsets=ethogram_base - count, colors=[colored_behaviors[behavior]], linelengths=linelengths)
            # Add a label for the behavior
            # ax_ethogram.text(fp_time[-1], ethogram_base - count, behavior, verticalalignment='center', color=colored_behaviors[behavior], fontsize=8)
            patch = mpatches.Patch(color=colored_behaviors[behavior], label=behavior)
            legend_patches.append(patch)
            count += 0.1
    
    # Adjust the ethogram plot's y-limits to tighten up the space around the plotted events
    ax_ethogram.set_ylim([ethogram_base - count, 0])
    ax_ethogram.set_yticks([])  # Hide y-axis ticks for the ethogram
    #set xticks
    ax_ethogram.set_xlabel('Time (s)')  # Add an x-axis label
    # Add legend to the ethogram axis
    ax_trace.legend(handles=legend_patches, loc='upper right', fontsize='small', frameon=False)
    
    # Hide the x-axis labels for the top plot (photometry trace)
    # ax_trace.tick_params(labelbottom=False, show = False)
    sub_folder = make_folder(f"{subject}", Representative_image)
    # Set up the layout so plots align nicely
    plt.tight_layout(pad = 0)
    plt.xlim(-fp_time.min(),fp_time.max())
    plt.show()
    # Save the figure to a file
    fig.savefig(f'{sub_folder}/{id_}_{trace_name}.svg')
    plt.close(fig)  # Close the figure to free up memory

# You would call `plot_ethogram_trace` within your code after processing your behavior data into `new_behav_dict` and preparing `fp_time`.
    
def plot_me(new_behav_dict, smooth_trace, fp_time, fp_start, fp_end, i):

    # color_list=['green', 'gold','red', 'lawngreen', 'dodgerblue',  'green', 'deeppink', 'aquamarine', 'deepskyblue','palevioletred', 'orange', 'mediumturquoise', 'red','sandybrown', 'orchid', 'springgreen', 'darkorchid', 'royalblue']
    color_list = ['dodgerblue', 'lawngreen', 'gold', 'magenta']
    # # create a color map with random colors
    colmap = matplotlib.colors.ListedColormap(color_list)


    #ETHOGRAM ALIGNED WITH NORMALIZED TRACE
    plt.figure(figsize = (10,4))

    # plt.align_ylabels()


    plt.plot(fp_time,smooth_trace, color = 'k')
    # etho_fig.label_outer()
    plt.ylim([smooth_trace.min()-(0.30*smooth_trace.min()), smooth_trace.max()+(0.30*smooth_trace.max())])
    plt.ylabel(r'$\Delta$F', fontsize = 20)
    plt.xlabel('Time (s)', fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    # plt.set_title(f'{i}')

    x_positions = np.arange(0,fp_time.max(),(100*fp_fps))
    # etho_fig.xticks(x_labels)  
    # print(new_behav_dict)
    plt.tight_layout()
    count = 0
    stagger = 0
    add_to_legend = []
    for key in new_behav_dict:
        add_to_legend.append(key)
        for k in range(len(new_behav_dict[key])):
            if new_behav_dict[key][k]==1:
                a=int(k)
                b=int(k+post_s*fp_fps)
                c=int(k-pre_s*fp_fps)
                plt.axvspan(fp_time[a], fp_time[b], color=color_list[count],  alpha=0.5, ymax = 0.4-stagger, label = key)
                # plt.axvspan(fp_time[a], fp_time[c], color=color_list[count], alpha=0.2, ymax = 0.2)
                plt.legend(fontsize = 20)
        count+=1
        stagger+=0.03
               

    plt.show()
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
        if m >= 0 and q+200 <= len(smooth_trace):
            baseline=[]
            event=[]
            base_start=[]
        #RETRIEVES VALUES FROM THE SMOOTH NORMALIZED TRACE
            for j, k in zip(smooth_trace[m:n], smooth_trace[i:q]): #Overshoot for AUC
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
                    # if j not in check_first:
                    placement.append(j)   
                if behav[o]==1 and behav[c]==1:
                    # print('yes')
                    j=floor((o/behav_fps)*fp_fps)
                    if j in l_base: # and j not in check_first:
                        placement.append(j)  
            else:
                if behav[o]==1 and behav[a]==0:
                    j=floor((o/behav_fps)*fp_fps)
                    placement.append(j) 
                elif behav[o]==0:
                    continue

    fix_place_2=placement[1:]; fix_place=placement[0:-2]

    for d,j in zip(fix_place, fix_place_2):
        if int(j) in range(d, (d+(pre_s*fp_fps))) or int(j-((pre_s-(floor(pre_s/2)))*fp_fps)) in range(d, (d+(floor(post_s/1.5)*fp_fps))):  #should prevent overlapping ROIs
            if j in placement:
                placement.remove(j)

    return placement

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
    if mini>0:
        behav_raw_=behav_raw[start-time_from_start*behav_fps: end+time_from_end*behav_fps]
        new_470=_470[int(fp_start-time_from_start*fp_fps):int(fp_end+time_from_end*fp_fps)]
        new_410=_410[int(fp_start-time_from_start*fp_fps):int(fp_end+time_from_end*fp_fps)]
    elif mini<0:
        # mini2=0-mini
        # mini3=0-(int(fp_start-time_from_start*fp_fps))
        behav_raw_=behav_raw[0: end+time_from_end*behav_fps]
        new_470=_470[0:int(fp_end+time_from_end*fp_fps)]
        new_410=_410[0:int(fp_end+time_from_end*fp_fps)]

    _470_=np.array([i for i in new_470])
    _410_=np.array([i for i in new_410])
    frame = np.array([i for i in range(len(_470_))])
    time = np.array([i/fp_fps for i in frame])

    behav_frames = [i for i in range(len(behav_raw_[behavior[0]]))]
    # TIME
    behav_time = np.array([i/behav_fps for i in behav_frames])

    return _470_,_410_, behav_raw_, behav_time, time, frame, behav_frames;

def takeoff(behav_time, fp_time, behavior, behav_fps, behav_raw):
    
    # fp_addon_offset=behav_duration-fp_time[-1]
    # #rationale: the behavior starts before photometry; but the behavior video decides when the photometry should start
    if behav_time.max()>fp_time.max():
        offset=(behav_time.max()-fp_time.max())
        cutit=(offset*behav_fps) ##Make sure it is rounding down
        # print(cutit)
        cutit=int(cutit)
        # print(cutit)
        # cutit=math.ceil(cut)
        # behav_cutit=int((offset*behav_fps)
        j=behav_raw[behavior[0]]
        behav_frames=np.array([i for i in range(len(j[cutit:]))])
        behav_times=[]
        for i in behav_frames:
            j=i/behav_fps
            behav_times.append(j)
        behav_times=np.array(behav_times)
        behav_raw=behav_raw[cutit:]

    return behav_times, behav_raw, cutit;

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

def photo_info(_470, _410, fp_frames):

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

    smooth_trace = smooth_trace_
    scaled_trace = smooth_trace*scale

    u = np.mean(smooth_trace)
    std = np.std(smooth_trace)
    zscore_trace = (smooth_trace-u)/std      
    
    plot_photo(fp_frames, _470, _410, fit1, fit3, smooth_trace)
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
        good_control = input("Remove isobesctic control? (y/n)")  
        if good_control.lower() == 'y':
            fp_times=[i/fp_fps for i in fp_frames]
            normalized=(_470-fit3)/fit3
            smooth_normal = smooth(_470, 59)  # must be an odd number
            smooth_trace_ = smooth_normal-smooth_normal.min()
            scale=1/smooth_trace_.max()
            smooth_trace = smooth_trace_
            scaled_trace = smooth_trace*scale

            u = np.mean(smooth_trace)
            std = np.std(smooth_trace)
            zscore_trace = (smooth_trace-u)/std 
            break
        if good_control.lower() == 'n' or good_control.lower() =='':
            fp_times=[i/fp_fps for i in fp_frames]
            normalized=(_470-fit3)/fit3
            smooth_normal = smooth(normalized, 59)  # must be an odd number
            smooth_trace_ = smooth_normal-smooth_normal.min()
            scale=1/smooth_trace_.max()
            smooth_trace = smooth_trace_
            scaled_trace = smooth_trace*scale

            u = np.mean(smooth_trace)
            std = np.std(smooth_trace)
            zscore_trace = (smooth_trace-u)/std 
            break

    # LIST OF PEAK POSITIONS
    x, _ = list(signal.find_peaks(smooth_trace, width=50, height=0.21))
    prom, l_base, r_base = list(signal.peak_prominences(smooth_trace, x, wlen=300))
    heights = smooth_trace[x]

    return smooth_trace, scaled_trace, zscore_trace, fp_frames, heights, l_base;

def open_files(fp_file, behav_file):

    try:

        fp_raw = pd.read_csv(fp_file)

        behav_raw = pd.read_excel(behav_file, header=[behavior_row-1], skiprows=[behavior_row])

        behavior, behav_frames, behav_time = behav_split(
            behav_raw, behav_fps)

        if lolo not in behavior:
            new_row=int(input("Behaviors are not in designated row. Try again: "))
            behav_raw = pd.read_excel(behav_file, header=[new_row-1], skiprows=[new_row])
                    # behavior_row=new_row
    except ValueError:
        print(f"Unable to open file(s)")
        exit()

    return fp_raw, behav_raw;

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

def plot_photo(fp_frames, _470, _410, fit1, fit3, smooth_trace):

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
    # ax[0].set_title(f'{project_id}')
    # ax[3].set(xlabel='Frame #', ylabel=r'F mean pixels')

    # 410 FIT TO 470
    ax[2].plot(fp_frames, _470)
    ax[2].plot(fp_frames, fit3, 'r')
    ax[2].set_title('410 fit over 470')
    ax[2].set(xlabel='Frame #', ylabel=r'F mean pixels')

    ax[3].plot(fp_frames, smooth_trace)
    plt.show(block=False)

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
