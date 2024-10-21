# AUTHOR: VICTORIA SEDWICK
# ADAPTED FROM NEUROfpMETRICS AND ILARIA CARTA (MATLAB)

# https://mat#plotlib.org/stable/gallery/sub#plots_axes_and_figures/sub#plots_demo.html
# https://stackoverflow.com/questions/11190735/python-mat#plotlib-superimpose-scatter-#plots

# <Where are your files located?>


fp_file = r"\\data.einsteinmed.edu\users\Autry Lab\Victoria Sedwick\2_cno(0)_meP_2024-06-25T20_55_21.CSV"
# NOTE <Load any ethovision file that has behavior names listed>
# example_behav_file = r"E:\Photometry-Fall2022\Pup Block\Virgin\Master file (Males)\Infanticidal Animals_trials\Pups\10\Raw data-MeP-CRFR2photo_pupblock1-Trial    31 (1).xlsx"
# <Frames per second for photo and ethovision data>
fp_fps = 20

#YES OR NOl
full_trace='yes'  #if no, the program will crop around the start and end times

#if full_trace = yes
crop_end=0; crop_front=0

time_start_at_drug_s = 10
time_start_in_drug_s  = int(23.5*60)


# Which roi/fiber signal do you want to #plot (0-2)?
z = 0
#PROCEED WITH CAUTION CHANGING ANYTHING AFTER THIS LINE
#MAKE SURE ALL PACKAGES ARE INSTALLED
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import exp
from pathlib import Path
import numpy as np
from scipy import stats, signal
from brokenaxes import brokenaxes


def main():
    global fp_file, fp_fps, z, time_start_at_drug_s, time_start_in_drug_s, crop_front, crop_end;

    fp_raw = pd.read_csv(fp_file)
    _470, _410 = fp_split(fp_raw, z, fp_fps)
    _470,_410, fp_time, fp_frames=cropendtime(_470, _410, crop_end) #cuts the last minute
    smooth_trace, scaled_trace, zscore_trace, fp_frames = photo_info(_470, _410, fp_frames)        
    trace_diary = {'Scaled_trace': scaled_trace, 'Zscore_trace': zscore_trace, 'Normalized trace': smooth_trace}
    for trace in trace_diary:
        segment1, segment2,time_s1, time_s2 = plot_first5_25(trace, trace_diary[trace], fp_time)
        plot_zscore_diff(segment1, segment2,time_s1, time_s2)
        plot_zExtraLabels(segment1, segment2,time_s1, time_s2)


# GET PHOTOMETRY VALUES FROM EXCEL SHEET
def fp_split(fp_raw, z, fp_fps):

    # Extract column headers from CSV to give user options
    roi = [i for i in fp_raw if 'Region' in i]
    # Indexing column headers, not values
    roi_ = roi[z]; led = fp_raw['LedState']; trace = fp_raw[roi_]

    c = 0
    # print(trace)
    _470_ = [j for i, j in zip(led, trace) if i == 6 and (c := c + 1) > (crop_front * fp_fps)]
    _410_ = [j for i,j in zip(led, trace) if i != 6 and (c := c + 1) > (crop_front * fp_fps)]
    
    if len(_470_) != len(_410_):
        if len(_470_) > len(_410_):
            _470_.remove(_470_[-1])
        elif len(_470_) < len(_410_):
            _410_.remove(_410_[-1])

    _470 = np.array(_470_)
    print(_470)
    _410 = np.array(_410_)

    frame = np.array([i for i in range(len(_470))])

    time = np.array([i/fp_fps for i in frame])

    return _470, _410

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

    return smooth_trace, scaled_trace, zscore_trace, fp_frames;


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

def cropendtime(_470, _410, crop_end):
    _470_=np.array([_470[i] for i in range(len(_470)) if i<(len(_470)-(crop_end*fp_fps))])
    _410_=np.array([_410[i] for i in range(len(_410)) if i<(len(_410)-(crop_end*fp_fps))])
    frame = np.array([i for i in range(len(_470_))])
    time = np.array([i/fp_fps for i in frame])

    return _470_,_410_,time, frame;

def plot_first5_25(trace_name, trace, fp_time):
    segment1 = trace[time_start_at_drug_s*fp_fps:(time_start_at_drug_s+300)*fp_fps]
    time_s1 = fp_time[time_start_at_drug_s*fp_fps:(time_start_at_drug_s+300)*fp_fps]

    segment2 = trace[time_start_in_drug_s*fp_fps:(time_start_in_drug_s+300)*fp_fps]
    time_s2 = fp_time[time_start_in_drug_s*fp_fps:(time_start_in_drug_s+300)*fp_fps]

    # Create subplots for the broken axis
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
    fig.subplots_adjust(wspace=0.1)  # Adjust the space between plots

    # Plot the first segment
    ax1.plot(time_s1, segment1, 'k', linewidth = 0.7)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Photometry Signal')
    ax1.set_title('First 5min')

    # Plot the second segment
    ax2.plot(time_s2, segment2, 'k', linewidth = 0.7)
    ax2.set_xlabel('Time (s)')
    ax2.set_title('20-25 minutes')  #make sure

    # Hide the spines between ax1 and ax2
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax2.yaxis.tick_right()

    # Add break lines
    d = .015  # size of the break marks
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the right axes
    ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    plt.suptitle('Photometry Trace with Broken Axis')
    # plt.show()

    return segment1, segment2, time_s1, time_s2

def plot_zscore_diff(segment1, segment2,time_s1, time_s2):
    u = np.mean(segment1)
    std = np.std(segment1)
    full_trace = np.concatenate((segment1, segment2))
    score_section = np.array([(i - u) / std for i in full_trace])

    fig = plt.figure(figsize = (9,5))
    bax = brokenaxes(xlims=((time_s1[0], time_s1[-1]), (time_s2[0], time_s2[-1])), hspace=0.05)

    # Plot the first segment
    bax.plot(time_s1, score_section[:len(time_s1)], 'k', linewidth=0.7)

    # Plot the second segment
    bax.plot(time_s2, score_section[len(time_s1):], 'k', linewidth=0.7)
    bax.set_ylabel('zF (relative to baseline)', size = 15)
    bax.set_xlabel('Time (s)', size = 15, labelpad = 25)
    bax.tick_params(axis='both', which='major', labelsize=12)

    bax.set_ylim([-2, 8])

    plt.show()



def plot_zExtraLabels(segment1, segment2,time_s1, time_s2):
    u = np.mean(segment1)
    std = np.std(segment1)
    full_trace = np.concatenate((segment1, segment2))
    score_section = np.array([(i - u) / std for i in full_trace])

    # Create subplots for the broken axis
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
    fig.subplots_adjust(wspace=0.1)  # Adjust the space between plots

    # Plot the first segment
    ax1.plot(time_s1, score_section[:len(time_s1)], 'k', linewidth = 0.7)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('zF (Relative to Baseline)')
    ax1.set_title('First 5min')

    # Plot the second segment
    ax2.plot(time_s2, score_section[len(time_s1):], 'k', linewidth = 0.7)
    ax2.set_xlabel('Time (s)')
    ax2.set_title('20-25 minutes')  #make sure

    # Hide the spines between ax1 and ax2
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax2.yaxis.set_visible(False)

    ax1.set_ylim([-2, 8])
    ax2.set_ylim([-2, 8])
    # Add break lines
    d = .015  # size of the break marks
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the right axes
    ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    # plt.suptitle('Photometry Trace with Broken Axis')
    plt.show()




if __name__ == "__main__":
    main()