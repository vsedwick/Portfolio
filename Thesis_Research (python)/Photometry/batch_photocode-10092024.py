
configuration_file = ''

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
import sys
import yaml
from pathlib import Path

def main(configuration_file):
    config = AnalysisParameters(load_config(configuration_file))

    ##NOTE Prompt for GUI to select behaviors to extract

    #Begin iterating through analysis folders
    analysis_folders = os.listdir(config.project_home)
    for subject_trial_id in analysis_folders:
        #Exclude processed value or video folders that may be in the project home directory
        exclude_folder = ['Behaviors', 'Videos', 'Summary']
        if not os.path.isdir(os.path.join(config.project_home, subject_trial_id)) or subject_trial_id in exclude_folder:
            continue
        else:
            #Identify and load raw files
            processor = Clean_and_ProcessFiles(config, subject_trial_id); 
            photo_470, photo_410, file_behavior_list, behav_raw_df = processor.process_files()


            #process, plot, and validate photometry plots
            #initialize new class
            #run the analysis function and return trace_diary
            trace_diary = #analysis function for photo specific class

            for i in trace_diary:

                processor.start_times; processor.end_times
                


       
            
class AnalysisParameters:
    def __init__(self, config):
        #Directory Settings
        self.project_home = Path(config['Directory_Information']['project_home']) 
        self.example_behavior_file = Path(config['Directory_Information']['example_behavior_file'])

        #Acquisition Settings
        self.behav_fps = config['Acquisition_Information']['behaviorvideo_fps']
        self.photo_fps = config['Acquisition_Information']['photometry_fps']
        self.pre_s = config['Acquisition_Information']['peri-baseline_seconds']
        self.post_s = config['Acquisition_Information']['peri-event_seconds']
        
        self.roi = config['Acquisition_Information']['roi']
        self.offset = config['Acquisition_Information']['offset']

        #Crop Settings        
        self.full_trace = config['Acquisition_Information']['full_trace']
        if self.full_trace.lower() == 'no':
            self.time_from_start = config['Acquisition_Information']['time_from_start_s']
            self.time_from_end = config['Acquisition_Information']['time_from_end_s']
        elif self.full_trace.lower()=='yes':
           self.crop_front = config['Acquisition_Information']['crop_front']
           self.crop_end = config['Acquisition_Information']['crop_end']
        else:
           raise ValueError("Parameter is missing from the configuration file: full_trace")

        #Behavior Parameters:
        self.behavior_row = config['Behavior_Parameters']['behavior_row']
        self.start_parameter = config['Behavior_Parameters']['start_parameter']
        self.min_duration = config['Behavior_Parameters']['minimum_accepted_duration_s']

        controls = config['Behavior_Parameters']['use_control_trial']
        if controls not in ['yes', 'no']:
            print("No control trial selected")
            self.controls = 'no'
        else:
            self.controls = controls
        
        self.perievent_limit = config['Behavior_Parameters']['limit_perievent_extraction']

        self.point_events = config['Behavior_Parameters']['point)events']

        compile_behavior = config['Behavior_Parameters']['Compile_behaviors']
        if compile_behavior.lower() not in ['yes', 'no']:
            self.compile_behavior = 'no'
        else:
            self.compile_behavior = compile_behavior
            if self.compile_behavior.lower() == 'yes':
                self.Groups = config['Behavior_Parameters']['Behavior_Groupings']
        
        add_zone = config['Behavior_Parameters']['Add_a_zone']
        if add_zone.lower() not in ['yes', 'no']:
            self.add_zone = 'no'
        else:
            self.add_zone = add_zone
            if self.add_zone.lower() == 'yes':
                self.zone = config['Behavior_Parameters']['Zone']
        
        self.report = []

class Clean_and_ProcessFiles:
    def __init__(self, config, subject_trial_id):
        self.config = config
        self.photo_raw = None
        self.behav_raw = None
        self.files_to_load = []
        self.photo_470 = []
        self.photo_410 = []
        self.subject_trial_id = subject_trial_id
        self.subject_trial_path = os.path.join(config.project_home, subject_trial_id)
        self.start_times = []
        self.end_times = []

    def open_files(self):
        self.files_to_load = os.listdir(self.subject_trial_path)
        print(f'Analyzing ID: {self.subject_trial_id}...'); print(f"{self.subject_trial_id}. Opening files")
        if len(self.files_to_load)>= 2:
            for i in self.files_to_load:
                try:
                    if i.lower().endswith('.csv'):
                        self.photo_raw = pd.read_csv(Path(os.path.join(self.subject_trial_path, i)))
                    if 'Raw data' in i and i.endswith('.xlsx') and '$' not in i:
                        self.behav_raw = pd.read_excel(Path(os.path.join(self.subject_trial_path, i)), header = [self.config.behavior_row-1], skiprows = [self.config.behavior_row])
                except (ValueError, FileNotFoundError, pd.errors.EmptyDataError):           
                    print(f"Unable to open file: {i}; skipping trial {self.subject_trial_id}")
                    self.config.report.append({self.subject_trial_id})
                    
                    continue
        else:
            self.config.report.append(self.subject_trial_id)
            print("Data not present. Skipping trial...")

        return self.photo_raw, self.behav_raw
    
    def split_photometry_values(self):

        photo_470 = []; photo_410 = []
        roi_index = [i for i in self.photo_raw if "Region" in i][self.config.roi]
        
        led_column = self.photo_raw['LedState']
        data_column = self.photo_raw[roi_index]
        counter = 0
        if self.config.use_full_trace == 'yes':
            for i, j in zip(led_column, data_column):
                if i==6 and counter > self.config.crop_front*self.config.photo_fps:
                    photo_470.append(j)
                elif i!= 6 and counter > self.config.crop_front*self.config.photo_fps:
                    photo_410.append(j)
                counter +=1
        else:
            photo_470 = [j for i, j in zip(led_column, data_column) if i == 6]
            photo_410 = [j for i, j in zip(led_column, data_column) if i != 6]
            
        
        if len(photo_470) != len(photo_410):
            if len(photo_470) > len(photo_410):
                photo_470 = photo_470[:len(photo_410)]
            else:
                photo_410 = photo_410[:len(photo_470)]

        self.photo_410 = np.array(photo_410); self.photo_470 = np.array(photo_470)

        return self.photo_410, self.photo_470


    def split_behavior_values(self):
        """
        Indexes the appropriate columns and rows to retreive information about the video data.
        
        Args:
            behaviors (list of str):
            behavior_frames (list of int):
            behavior_time (array):
        Return:
            behaviors
            ....
            """
        #get behavior names; based on Noldus output
        self.behaviors = [i for i in self.behav_raw][7:-1]

        return self.behaviors
   
    def validate_files(self):
        if self.config.start_parameter not in self.behaviors:
            self.config.report.append(self.subject_trial_id)
            return False
        else:
            return True
    
    def set_times_equal(self): #might have 1 sec difference

        behav_time = get_time_array(self.behav_raw, self.config.behav_fps)
        photo_time = get_time_array(self.photo_470, self.config.photo_fps)

        if behav_time.max() + 1 > photo_time.max():
            cutit = floor(self.config.offset)

            j = self.behav_raw[self.behaviors[0]]

            self.behav_raw = self.behav_raw[cutit:]
        
        return self.behav_raw

    def determine_start_times(self):
        start_time_array = np.array(self.behav_raw[self.config.start_parameter])
        start_time_placements = [int(i) for i in range(len(start_time_array)) if start_time_array[i] == 1]
        try:
            if self.config.use_controls == 'yes':
                if len(start_time_placements) < 4:
                    self.start_times = start_time_placements[0]; self.end_times = start_time_placements [1]
                    self.start_times = start_time_placements[2]; self.end_times = start_time_placements [3]
                else:
                    self.start_times = start_time_placements[2]; self.end_times = start_time_placements [3]
            elif self.config.use_controls == 'no': self.start_times = start_time_placements[0]; self.end_times = start_time_placements [1]
        except IndexError:
            print("The start parameter is not indicated for this subject")
            self.config.report.append(self.subject_trial_id)
            return

    def crop_file_ends(self):
        if self.config.full_trace == 'yes':
            raise ValueError("Check the 'full_trace' paramter in the configuration file. Make sure it says 'no'.") 
        else:
            crop_index = int(self.config.crop_end * self.config.photo_fps)

            self.photo_470 = self.photo_470[:-crop_index]
            self.photo_410 = self.photo_410[:-crop_index]
            self.behav_raw = self.behav_raw[:-int(self.config.crop_end * self.config.behav_fps)]

    def crop_around_start_param(self):
        if self.config.full_trace == 'yes':
            photo_start = int(self.start_times[0]/self.config.behav_fps*self.config.photo_fps)
            photo_end = int(self.end_times[-1]/self.config.behav_fps*self.config.photo_fps)

            pre_window_start = max(0, self.start_time[0] - self.config.time_from_start * self.behav_fps)
            post_end_limit = min(len(self.photo_470), photo_end + self.config.time_from_end * self.config.photo_fps)

            self.behav_raw = self.behav_raw[pre_window_start : self.end_times[-1] + self.config.time_from_end * self.config.behav_fps]
            self.photo_470 = self.photo_470[max(0, photo_start - self.config.time_from_start * self.behav_fps) : post_end_limit]
            self.photo_410 = self.photo_410[max(0, photo_start - self.config.time_from_start * self.behav_fps) : post_end_limit]

        else:
            raise ValueError("Check the 'full_trace' paramter in the configuration file. Make sure it says 'yes'.") 
    
    def correct_led_assignment(self):
        # Normalize both signals
        check_photo_410 = (self.photo_410[1000:] - self.photo_410[1000:].min()) / (self.photo_410[1000:].max() - self.photo_410[1000:].min())
        check_photo_470 = (self.photo_470[1000:] - self.photo_470[1000:].min()) / (self.photo_470[1000:].max() - self.photo_470[1000:].min())

        # Compute the range or variance of both signals
        range_410 = check_photo_410.max() - check_photo_470.min()
        range_470 = check_photo_470.max() - check_photo_410.min()

        # Detect if signal swap is needed based on range (isosbestic signal should have a smaller range)
        if range_410 > range_470:
            # Swap signals if needed
            self.photo_410, self.photo_470 = self.photo_470, self.photo_410

    def process_files(self):
            self.open_files()
            self.split_behavior_values()

            if not self.validate_files():
                return False
            else:
                self.split_photometry_values()

                self.correct_led_assignment()

                behav_time, photo_time = get_time_array(self.behav_raw, self.config.behav_fps), get_time_array(self.photo_470, self.config.photo_fps)

                #Correct frame rate
                print(f"Correcting frames per second. Current fps: {float(len(self.photo_470)/behav_time.max())}")
                self.correct_fps()

                

                print(f"Corrected frames per second: {len(self.photo_470)/max(photo_time)}")
                
                #prepare to crop
                print(f"Time difference between video and trace before cropping: {behav_time.max()}, {photo_time.max()}")
                if self.config.full_trace == 'yes':
                    self.crop_file_ends()
                else:
                    self.crop_around_start_param()

                new_photo_time = get_time_array(self.photo_470, self.config.photo_fps)
                behav_time = get_time_array(self.behav_raw, self.config.behav_fps)
                
                print(f"Time difference between video and trace after cropping: {new_behav_time.max()}, {new_photo_time.max()}")

                return self.photo_410, self.photo_410, self.behaviors, self.behav_raw
    
    def correct_fps(self): #NOTE fix
        """
        Matches the photometry frames per second to the video for consistency purposes
        """

        #Find current frames per second for photometry data
        input_fps = int(len(self.photo_470)/self.behavior_time)
        #Identify how many frames need to be added per second
        total_frames_togenerate = int(self.behavior_time*self.config.photo_fps)
        #interpolate values
        interp_index = np.linspace(0, len(self.photo_470) - 1, total_frames_togenerate)
        self.photo_470 = np.interp(interp_index, np.arange(len(self.photo_470)), self.photo_470)
        self. photo_410 = np.interp(interp_index, np.arange(len(self.photo_410)), self.photo_410)

class SmoothPlotnValidate_Photo:
    def __init__(self, config, photo_470, photo_410, trial_id):
        self.photo_470 = photo_470
        self.photo_410 = photo_410
        self.config = config
        self.id = trial_id
        self.fit1 = []
        self.fit3 = []
        self.remove_isocesctic = False
        self.heights = []
        self.trace_diary = None
        self.photo_frames = []
        self.proper_assignment = True
    def exp2(self, x, a, b, c, d): # NOTE check what these variables are
        return a*exp(b*x) + c*exp(d*x)
    
    def curve_fitting (self, photo_470, photo_410):
        self.photo_frames = get_frame_array(photo_470)
        
        popt, pcov = curve_fit(self.exp2, self.photo_frames, photo_410, maxfev = 500000, p0 = [0.01, 1e-7, 1e-5, 1e-9])
        self.fit1 = self.exp2(self.self.photo_frames, *popt)

        A = np.vstack([self.fit1, np.ones(len(self.fit1))]).T # NOTE rename variable
        slope = np.linalg.lstsq(A, photo_470, rcond = None)[0]
        fit2 = self.fit1 * slope[0] + slope[1]

        # FIT 410 OVER 470
        fit3 = stats.linregress(photo_410, photo_470)
        self.fit3 = fit3.intercept+fit3.slope * photo_410

        return self.fit1, fit2, self.fit3

    # SMOOTH AND NORMALIZE PHOTOMETRY TRACE
    def smooth(a, WSZ):  # 
        """
        Formula to smooth wave forms obtained from https://stackoverflow.com/questions/40443020/matlabs-smooth-implementation-n-point-moving-average-in-numpy-python.
        
        Args:
            a (array): NumPy 1-D array containing the data to be smoothed (e.g. photo_470).
            WSZ (int??): smoothing window size; should be an odd number
        
        Returns:
            ???
            """

        out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid')/WSZ
        r = np.arange(1, WSZ-1, 2)
        start = np.cumsum(a[:WSZ-1])[::2]/r
        stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
        return np.concatenate((start, out0, stop))

    def normalize_and_zscore_trace(self, photo_470):
        
        fitted_trace = (photo_470-self.fit3)/self.fit3
        if self.remove_isobesctic:
            smoothed_and_fitted_trace = self.smooth(photo_470, 59)
        else:
            smoothed_and_fitted_trace = self.smooth(fitted_trace, 59)
        self.normalized_trace = smoothed_and_fitted_trace - smoothed_and_fitted_trace.min()
        
        u = np.mean(self.normalized_trace)
        std = np.std(self.normalized_trace)
        self.zF_trace = (self.normalized_trace-u)/std

        self.scaled_normalized_trace = self.normalized_trace * (1/self.normalized_trace.max())

        self.trace_diary = {'Normalized_trace': self.normalized_trace, 'Scaled_trace': self.scaled_normalized_trace, 'Zscore_trace': self.zF_trace}

        return self.trace_diary

        ###NOTE consider adding a validation prompt.
    
    def identify_peak_positions(self, trace):
        # LIST OF PEAK POSITIONS
        x, _ = list(signal.find_peaks(trace, width=50, height=0.21))
        prom, l_base, r_base = list(signal.peak_prominences(trace, x, wlen=300))
        self.heights = trace[x]
        return self.heights
    
    def plot_photo_fittings(self, trace_in_diary):
        photo_frames = get_frame_array(self.photo_470)

        fig, ax = plt.subplots(4)
        #plot 1: raw trace
        ax[0].plot(photo_frames, self.photo_470)
        ax[0].set_title('Raw 47\0')
        ax[0].set(xlabel='Frame #', ylabel=r'$\Delta$F/F')

        # BIEXPONENTIAL FIT
        ax[1].set_ylim([self.photo_410.min(), self.photo_410.max()])
        ax[1].plot(photo_frames, self.photo_410)
        ax[1].plot(photo_frames, self.fit1, 'r')
        # rationale: https://onlinelibrary.wiley.com/doi/10.1002/mma.7396
        ax[1].set_title('410 with biexponential fit')
        ax[1].set(xlabel='Frame #', ylabel=r'$\Delta$F/F')
        ax[0].set_title(f'{self.trial_id}')
        # ax[3].set(xlabel='Frame #', ylabel=r'F mean pixels')

        # 410 FIT TO 470
        ax[2].plot(photo_frames, self.photo_470)
        ax[2].plot(photo_frames, self.fit3, 'r')
        ax[2].set_title('410 fit over 470')
        ax[2].set(xlabel='Frame #', ylabel=r'F mean pixels')

        ax[3].plot(photo_frames, trace_in_diary)
        plt.show(block=False)
        

    def get_clean_traces(self, photo_470, photo_410):
        
        self.curve_fitting(photo_470, photo_410)

        trace_diary = self.normalize_and_zscore_trace(photo_470)
        #plot
        self.plot_photo_fittings(trace_diary[0])

        #In case previous led correction was not successful
        while True:
            proper_led_assignment = input("Are 470 and 410 LEDs properly assigned? ('y'/ENTER or 'n'; input 's' to skip trial) \n NOTE: can silence this prompt at line _____: ")
            plt.close()
            if proper_led_assignment.lower() in ['n', 'no']:
                self.photo_410, self.photo_470 = photo_470, photo_410
                self.curve_fitting(photo_470, photo_410)
                trace_diary = self.normalize_and_zscore_trace(photo_470)
                self.plot_photo_fittings(trace_diary[0])
                break

            if proper_led_assignment.lower() == 's':
                self.config.report.append(self.trial_id)
                print(f"Skipping {self.trial_id}")
                return

        #In case isobesctic signal is not good.
        while True:
            good_isobesctic = input("Do you need to remove the isobesctic control signal? ('y'/'n'; press ENTER to skip): ")
            plt.close()
            if good_isobesctic.lower() in ['y', 'yes']: 
                self.remove_isocesctic = True
                trace_diary = self.normalize_and_zscore_trace(self.photo_470)
                break
            else:
                break
        
        #validate
        #get peak positions
        self.identify_peak_positions(trace_diary[0])

        return trace_diary

class Restructure_Behavior_Data:
    def __init__(self, behav_raw, events_to_extract, use_zones):
        """
        Description #NOTE
        
        Attributes:
        events_to_extract (lst of str): list of scored event information you wish to extract from photometry trace.
        event_dictionary (dict): Newly restructed event data that allows ethogram plotting
        """
        self.events_to_extract = events_to_extract
        self.behav_raw = behav_raw
        self.event_dictionary = None
        self.use_zones = use_zones

    def create_dictionary(self):
        for event in self.events_to_extract:
            bool_array = np.array(self.behav_raw[event], dtype = bool)
            self.event_dictionary.update({f"{event}": bool_array})

        return self.event_dictionary
    
    def define_zoned_behaviors(self, event_dictionary):
        for event in self.events_to_extract:
            in_zone = []; out_of_zone = []
            for event_value, zone_value in zip(self.behav_raw[event_value], self.behav_raw[zone_value]):
                if event_value == 1 and zone_value == 1:
                    out_of_zone.append(0)
                    in_zone.append(1)
                if event_value == 1 and zone_value == 0:
                    out_of_zone.append(1)
                    in_zone.append(0)
                if (event_value == 0 and zone_value == 0) or (event_value == 0 and zone_value == 1):
                    out_of_zone.append(0)
                    in_zone.append(0)                    
            in_zone_array = np.array(in_zone, dtype = bool)
            out_zone_array = np.array(out_of_zone, dtype = bool)

            in_zone_label = f"{event}_InZone"
            out_zone_label = f"{event}_OutZone" 

            self.events_to_extract.append([in_zone_label]); self.events_to_extract.append([out_zone_label])

            event_dictionary.update({f"{in_zone_label}": in_zone_array,
                                     f"{out_zone_label}": out_zone_array})
            self.event_dictionary = event_dictionary
            
            return self.event_dictonary
        
        def restructure(self):
            self.create_dictionary()
            if self.use_zones in ['y']


        


def load_config(variables_config):
    """ 
    Loads the configuration file. can be specified as an argument or defined in the script.
    
    Args: 
        variables_config (str): The file path that directs to the config file. 
        
    Returns:
        config (dict): A dictionary of the configuration file that contains the analysis parameters.

    """
    with open(variables_config, 'r') as file:
        config = yaml.safe_load(file)
    return config

def get_frame_array(array1):
    return [i for i in range(len(array1))]

def get_time_array(array2, fps):
    return [i/fps for i in range(len(array2))]

def transform_behavior_dataframe(behav_raw, events_to_extract, config): #NOTE behaviors to score are not yet defined.
    for i in events_to_extract:
        event_frame_occurance = [int(frame/config.behav_fps/config.photo_fps) for frame, value in zip(range(len(behav_raw[i])), behav_raw[i]) if value == 1]
        


if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(configuration_file)


        

