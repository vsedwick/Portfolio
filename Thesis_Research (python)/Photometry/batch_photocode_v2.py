
configuration_file = r'C:\Users\sedwi\Desktop\Portfolio\Thesis_Research (python)\Photometry\config.yaml'

# NOTE PROCEED WITH CAUTION CHANGING ANYTHING AFTER THIS LINE
# MAKE SURE ALL PACKAGES ARE INSTALLED WITH 'pip install [package]'
import os
import yaml
from math import floor
from numpy import mean
from scipy import stats, signal
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path
import yaml
from utils.logger import logger



            
class AnalysisParameters:
    """
    Loads analysis parameters from the YAML configuraton file.
    Refer to the YAML file for detailed descriptions of specific attributes.

    Attributes:
        project_home (Path): Path to the project's home directory.

        silent_mode (bool): Whether to run the acquisition in silent mode.
        full_trace (bool): Whether to use the full photometry trace or crop it.
        controls (bool): Whether to use control trials.
        compile_behavior (bool): Whether to compile behaviors.
        add_zone (bool): Whether to include a custom zone for analysis.

        behav_fps (float): Frames per second for the behavior video.
        photo_fps (float): Frames per second for the photometry data.
        pre_s (int): Time before an event for peri-event baseline in seconds.
        post_s (int): Time after an event for peri-event baseline in seconds.
        crop_front (int): Front crop value if full trace is used.
        crop_end (int): End crop value if full trace is used.
        time_from_start (int): Time in seconds from the start to crop, if full trace is not used.
        time_from_end (int): Time in seconds from the end to crop, if full trace is not used.
        min_duration (float): Minimum accepted duration for a behavior event in seconds.

        roi (int): LED region of interest for the analysis.
        point_events (list): List of point events for behavior analysis.
        Groups (dict): Groupings of behaviors for analysis if compile_behavior is True.

        report (list): List that gets updated throughout the program for reporting.
        trial_id (int): Identifier for the current trial.
        event_name (str): Name of the current event.
        stage (str): Current stage of the analysis process.

     Methods:
        _validate_path(path, param_name): Validates the existence of a path.
        _validate_boolean(value, param_name): Validates that a value is boolean.
        _clean_dict(d): Cleans a dictionary by removing empty or None values.
"""

    #@logger("log.txt")
    def __init__(self, config):
        """
        Initializes the variables from the provided YAML file.
        
        Args:
            config (dict): Configuration dictionary.

        """
        #Directory Settings
        
        self.project_home = self._validate_path(Path(config['Directory_Information']['project_home']), 'project_home')
        #Acquisition Settings
        self.silent_mode = config['Acquisition_Information']['silent_mode']

        self.behav_fps = config['Acquisition_Information']['behaviorvideo_fps']
        self.photo_fps = config['Acquisition_Information']['photometry_fps']
        self.pre_s = config['Acquisition_Information']['peri-baseline_seconds']
        self.post_s = config['Acquisition_Information']['peri-event_seconds']
        
        self.roi = config['Acquisition_Information']['roi']
        self.offset = config['Acquisition_Information']['offset']

        #Crop Settings        
        self.full_trace = self._validate_boolean(config['Acquisition_Information']['full_trace'], 'Acquisition_Information: full_trace' )
        if not self.full_trace:
            self.time_from_start = config['Acquisition_Information']['time_from_start_s']
            self.time_from_end = config['Acquisition_Information']['time_from_end_s']
        elif self.full_trace:
           self.crop_front = config['Acquisition_Information']['crop_front']
           self.crop_end = config['Acquisition_Information']['crop_end']

        #Behavior Parameters:
        self.behavior_row = config['Behavior_Parameters']['behavior_row']
        self.start_parameter = config['Behavior_Parameters']['start_parameter']
        self.min_duration = config['Behavior_Parameters']['minimum_accepted_duration_s']

        self.controls = self._validate_boolean(config['Behavior_Parameters']['use_control_trial'], 'Acquisition_Information: full_trace')
        
        self.perievent_limit = config['Behavior_Parameters']['limit_perievent_extraction']

        self.point_events = config['Behavior_Parameters']['point_events']

        #Events to score

        self.events_to_score = [i for i in config['Behavior_Parameters']['Behaviors_to_Score'] if i is not None]

        print("Events to Score: ", self.events_to_score)
        self.compile_behavior = self._validate_boolean(config['Behavior_Parameters']['Compile_behaviors'], 'Behavior_Parameters: Compile_behaviors')
        if not self.compile_behavior:
            self.Groups = None
        else:
            self.Groups = self._clean_dict(config['Behavior_Parameters']['Behavior_Groupings'])
        print(self.Groups)
        self.add_zone = self._validate_boolean(config['Behavior_Parameters']['Add_a_zone'], 'Behavior_Parameters, Add_a_zone')
        if not self.add_zone:
            self.zone = None
        else:
            self.zone = config['Behavior_Parameters']['Zone']

        #Updated throughout the program
        self.report = []
        self.trial_id = None
        self.event_name = None
        self.subject_path = None
        self.stage = None

    #@logger("log.txt")
    def _validate_path(self, path, param_name):
        """
        Validates the existence of the provided directory path.

        Args:
            path (Path): Path to validate.
            param_name (str): Name of the parameter being validated.

        Returns:
            Path: Validated path if it exists.
        
        Raises:
            FileNotFoundError: If the path does not exist.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{param_name} directory '{path}' does not exist.")
        else:
            return path

    #@logger("log.txt")
    def _validate_boolean(self, value, param_name):
        """
        Validates that a parameter is a boolean.

        Args:
            value (any): Value to validate.
            param_name (str): Parameter name for error reporting.

        Returns:
            bool: True if the value is a boolean.
        
        Raises:
            ValueError: If the value is not a boolean.
        """

        if value is None:
            return False
        if not isinstance(value, bool):
            raise ValueError(f"The parameter '{param_name}' must be a boolean.")
        return value
    
    def _clean_dict(self, d):
        """
        Cleans up a dictionary by removing empty or None values.

        Args:
            d (dict): Dictionary to clean.

        Returns:
            dict: Cleaned dictionary with non-empty values.
        """

        keys_to_remove = []  # List of keys to remove
        
        for key, value in d.items():
            if isinstance(value, list):
                # Remove None values from lists
                d[key] = [item for item in value if item is not None]
                # Mark for removal if the list is now empty
                if not d[key]:  
                    keys_to_remove.append(key)
            
            elif isinstance(value, dict):
                # Recursively clean nested dictionaries
                d[key] = clean_dict(value)
                # Mark for removal if the dictionary is now empty
                if not d[key]:  
                    keys_to_remove.append(key)
            
            elif value is None:
                # Mark for removal if the value is None
                keys_to_remove.append(key)

        # Remove keys with empty values
        for key in keys_to_remove:
            del d[key]

        return d
class Clean_and_ProcessFiles:
    #@logger("log.txt")
    def __init__(self, config):
        self.config = config
        self.photo_raw = None
        self.behav_raw = None
        self.files_to_load = []
        self.photo_470 = []
        self.photo_410 = []
        self.subject_trial_id = self.config.trial_id
        self.subject_trial_path = os.path.join(config.project_home, self.subject_trial_id)
        self.start_times = []
        self.end_times = []
        self.behav_time = []
        self.photo_time = []
        self.behaviors = []

    #@logger("log.txt")
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
                        check_header = list(self.behav_raw.head())
                        #specific to Noldus
                        if not check_header[0] == 'Trial time':
                            print("The provided behavior row is incorrect, finding the correct one...")
                            self.behav_raw = self.find_header_row(Path(os.path.join(self.subject_trial_path, i)))

                except (ValueError, FileNotFoundError, pd.errors.EmptyDataError):           
                    print(f"Unable to open file: {i}; skipping trial {self.subject_trial_id}")
                    self.config.report.append({self.subject_trial_id})
                    
                    return
        else:
            self.config.report.append(self.subject_trial_id)
            print("Data not present. Skipping trial...")

        return self.photo_raw, self.behav_raw
    
    #@logger("log.txt")
    def find_header_row(self, file):
        df = pd.read_excel(file, header = None)
        
        #specific to Noldus
        search_keywords = ['Trial time', 'Recording time', 'X center', 'Y center', 'Area', 'Areachange']
        proper_header = df.apply(lambda row: row.astype(str).str.contains('|'.join(search_keywords)).any(), axis=1)

        header_index = df[proper_header].index[0]
        self.config.behavior_row = header_index+1

        return pd.read_excel(file, header = [header_index], skiprows = [header_index + 1])
    
    #@logger("log.txt")
    def split_photometry_values(self):

        photo_470 = []; photo_410 = []
        roi_index = [i for i in self.photo_raw if "Region" in i][self.config.roi]
        
        led_column = self.photo_raw['LedState']
        data_column = self.photo_raw[roi_index]
        counter = 0
        if self.config.full_trace:
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

        self.photo_410 = np.array(photo_410) 
        self.photo_470 = np.array(photo_470)

        return self.photo_410, self.photo_470

    #@logger("log.txt")
    def split_behavior_values(self):
        """
        Indexes the appropriate columns and rows to retreive information about the video data.
        
        Args:
            behaviors (list of str):
        Return:
            behaviors
            ....
            """
        #get behavior names; based on Noldus output
        self.behaviors = [i for i in self.behav_raw][7:-1]

        return self.behaviors
   
    #@logger("log.txt")
    def validate_files(self):
        if self.config.start_parameter not in self.behaviors:
            self.config.report.append(self.subject_trial_id)
            return False
        else:
            return True
    
    #@logger("log.txt")
    def trim_beginning(self): #might have 1 sec difference

        behav_time = get_time_array(self.behav_raw, self.config.behav_fps)
        photo_time = get_time_array(self.photo_470, self.config.photo_fps)

        if behav_time.max() + 1 > photo_time.max():
            diff = behav_time.max() - photo_time.max()
            cutit = floor(diff * self.config.behav_fps)

            self.behav_raw = self.behav_raw[cutit:]
        
        return self.behav_raw

    #@logger("log.txt")
    def determine_start_times(self, behav_raw):
        self.start_times = []
        self.end_times = []
        try:
            start_time_array = np.array(behav_raw[self.config.start_parameter])
            start_time_placements = [int(i) for i in range(len(start_time_array)) if start_time_array[i] == 1]
        except ValueError:
            print("The start parameter is not indicated for this subject")
            self.config.report.append(self.subject_trial_id)
        # try:
        if self.config.controls:
            if len(start_time_placements) == 4:
                self.start_times.append(start_time_placements[0]); self.end_times.append(start_time_placements[1])
                self.start_times.append(start_time_placements[2]); self.end_times.append(start_time_placements[3])
            else:
                raise KeyError("'use_control_trial' is set to 'on'. There is not enough 'start' markers.")
        elif not self.config.controls: 
            self.start_times.append(start_time_placements[0]); self.end_times.append(start_time_placements[-1])

        return self.start_times, self.end_times

    #@logger("log.txt")
    def crop_file_ends(self):
        if not self.config.full_trace:
            raise ValueError("Check the 'full_trace' paramter in the configuration file. Make sure it says 'no'.") 
        else:
            crop_index = int(self.config.crop_end * self.config.photo_fps)
 
            self.photo_470 = self.photo_470[:-crop_index]
            self.photo_410 = self.photo_410[:-crop_index]
            self.behav_raw = self.behav_raw[:-int(self.config.crop_end * self.config.behav_fps)]

    #@logger("log.txt")
    def crop_around_start_param(self):
        if not self.config.full_trace:
            # Converting start and end times to photometry frames
            start = int(self.start_times[0])
            end = int(self.end_times[-1])

            # Convert to photometry frames
            photo_start = int(start / self.config.behav_fps * self.config.photo_fps)
            photo_end = int(end / self.config.behav_fps * self.config.photo_fps)

            # Calculate the check for boundaries
            check = start - self.config.time_from_start * self.config.behav_fps

            # Ensure cropping does not exceed array sizes
            if check > 0:
                pre_window_start = max(0, start - self.config.time_from_start * self.config.behav_fps)
                post_end_limit = min(len(self.photo_470), photo_end + int(self.config.time_from_end * self.config.photo_fps))

                # Adjust crop limits for behavior data
                behav_end = min(len(self.behav_raw), end + int(self.config.time_from_end * self.config.behav_fps))
                self.behav_raw = self.behav_raw[pre_window_start:behav_end]

                # Crop photometry data with boundary checks
                photo_start_cropped = max(0, photo_start - int(self.config.time_from_start * self.config.photo_fps))
                self.photo_470 = self.photo_470[photo_start_cropped:post_end_limit]
                self.photo_410 = self.photo_410[photo_start_cropped:post_end_limit]

            elif check < 0:
                # Handle case where start time is too early to subtract window
                behav_end = min(len(self.behav_raw), end + int(self.config.time_from_end * self.config.behav_fps))
                self.behav_raw = self.behav_raw[:behav_end]

                photo_end_limit = min(len(self.photo_470), photo_end + int(self.config.time_from_end * self.config.photo_fps))
                self.photo_470 = self.photo_470[:photo_end_limit]
                self.photo_410 = self.photo_410[:photo_end_limit]

            elif len(self.photo_470) < (photo_end + int(self.config.time_from_end * self.config.photo_fps)):
                # Handle case where photo trace is shorter than expected
                photo_start_cropped = max(0, photo_start - int(self.config.time_from_start * self.config.photo_fps))
                self.behav_raw = self.behav_raw[start - self.config.time_from_start * self.config.behav_fps:]
                self.photo_470 = self.photo_470[photo_start_cropped:]
                self.photo_410 = self.photo_410[photo_start_cropped:]

            else:
                raise ValueError("Check the 'full_trace' parameter in the configuration file. Make sure it says 'yes'.")


    #@logger("log.txt")
    def correct_led_assignment(self):
        # Normalize both signals
        check_photo_410 = (self.photo_410 - self.photo_410.min()) / (self.photo_410.max() - self.photo_410.min())
        check_photo_470 = (self.photo_470 - self.photo_470.min()) / (self.photo_470.max() - self.photo_470.min())

        # Calculate the signal-to-noise ratio (SNR) for both signals
        snr_410 = np.mean(check_photo_410) / np.std(check_photo_410)
        snr_470 = np.mean(check_photo_470) / np.std(check_photo_470)
        snr_diff = snr_410 / snr_470
        print(snr_diff)
        # Define thresholds for a significant difference
        snr_threshold = 1.5  # Adjust based on typical SNR range
        variance_threshold = 1.2  # Adjust based on variance ranges

        # Check if the 470 signal has a significantly higher SNR (as expected for the clear signal)
        if snr_410 > snr_470 and snr_diff > snr_threshold:
            # Swap signals if the 410 signal appears clearer than expected
            print(f"switching LEDs, 410: {snr_410}, 470: {snr_470}")
            self.photo_410, self.photo_470 = self.photo_470, self.photo_410
    # As a final check, ensure the 470 signal shows more variation (indicative of a clear signal)
        variance_410 = np.var(check_photo_410)
        variance_470 = np.var(check_photo_470)
        variance_diff = variance_470 / variance_410
        print("variance diff", variance_diff)

        if variance_410 > variance_470 and variance_diff > variance_threshold:
            print("switching LEDs")
            # Swap based on variance if the 410 signal is showing unexpected variation
            self.photo_410, self.photo_470 = self.photo_470, self.photo_410
    
    #@logger("log.txt")
    def correct_fps(self): 
        """
        Matches the photometry frames per second to the video for consistency purposes
        """

        #Find current frames per second for photometry data
        input_fps = float(len(self.photo_470)/self.behav_time.max())
        # Only correct if the rounded input FPS doesn't match the expected FPS
        if round(input_fps, 1) != self.config.photo_fps:
            print(f"Correcting frames per second. Current fps: {input_fps}")
            
            total_frames_togenerate = int(self.behav_time.max() * self.config.photo_fps)
            
            # Interpolate values
            interp_index = np.linspace(0, len(self.photo_470) - 1, total_frames_togenerate)
            self.photo_470 = np.interp(interp_index, np.arange(len(self.photo_470)), self.photo_470)
            self.photo_410 = np.interp(interp_index, np.arange(len(self.photo_410)), self.photo_410)

            output_fps = float(len(self.photo_470)/self.behav_time.max())
            
            print(f"Corrected frames per second: {output_fps}")

    #@logger("log.txt")
    def process_files(self):
            self.open_files()
            self.split_behavior_values()

            if not self.validate_files():
                return False
            else:
                self.split_photometry_values()

                self.correct_led_assignment()
                self.trim_beginning()

                self.determine_start_times(self.behav_raw)

                self.behav_time, self.photo_time = get_time_array(self.behav_raw, self.config.behav_fps), get_time_array(self.photo_470, self.config.photo_fps)

                #prepare to crop
                print(f"Time difference between video and trace before cropping: {self.behav_time.max()}, {self.photo_time.max()}")
                if self.config.full_trace:
                    self.crop_file_ends()
                else:
                    print("cropping around start param")
                    self.crop_around_start_param()
                
                self.determine_start_times(self.behav_raw)

                self.photo_time = get_time_array(self.photo_470, self.config.photo_fps)
                self.behav_time = get_time_array(self.behav_raw, self.config.behav_fps)
                # sys.exit()
                
                print(f"Time difference between video and trace after cropping: {self.behav_time.max()}, {self.photo_time.max()}")

                #Correct frame rate
                self.correct_fps()

                return self.photo_470, self.photo_410, self.behaviors, self.behav_raw
    

class SmoothPlotnValidate_Photo:
    def __init__(self, config, photo_470, photo_410):
        self.photo_470 = photo_470
        self.photo_410 = photo_410
        self.config = config
        self.fit1 = []
        self.fit3 = []
        self.remove_isosbestic = False
        self.trace_diary = None
        self.photo_frames = []
        self.proper_assignment = True
        self.trial_id = config.trial_id

    
    def fit_that_curve (self, photo_470, photo_410):
        self.photo_frames = get_frame_array(photo_470)
        
        popt, pcov = curve_fit(exp2, self.photo_frames, photo_410, maxfev = 500000, p0 = [0.01, 1e-7, 1e-5, 1e-9])
        self.fit1 = exp2(self.photo_frames, *popt)

        A = np.vstack([self.fit1, np.ones(len(self.fit1))]).T # NOTE rename variable

        slope = np.linalg.lstsq(A, self.photo_470, rcond = None)[0]
        fit2 = self.fit1 * slope[0] + slope[1]

        # FIT 410 OVER 470
        fit3_ = stats.linregress(self.photo_410, self.photo_470)
        self.fit3 = fit3_.intercept + fit3_.slope * self.photo_410
        return self.fit1, fit2, self.fit3

    # SMOOTH AND NORMALIZE PHOTOMETRY TRACE
    def smooth(self, a, WSZ):  # 
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

    #@logger("log.txt")
    def normalize_and_zscore_trace(self, photo_470):

        fitted_trace = np.array((photo_470-self.fit3)/self.fit3)

        if self.remove_isosbestic or fitted_trace.all() == 0:
            print("Removing isosbestic signal")
            smoothed_and_fitted_trace = self.smooth(photo_470, 59)
        else:
            smoothed_and_fitted_trace = self.smooth(fitted_trace, 59)
        self.normalized_trace = smoothed_and_fitted_trace - smoothed_and_fitted_trace.min()
        
        u = np.mean(self.normalized_trace)
        std = np.std(self.normalized_trace)
        
        self.zF_trace = (self.normalized_trace-u)/std

        self.scaled_normalized_trace = self.normalized_trace * (1/self.normalized_trace.max())
        self.trace_diary = {'Normalized_trace': self.normalized_trace, 
                             'Scaled_trace': self.scaled_normalized_trace, 
                             'Zscore_trace': self.zF_trace}

        return self.trace_diary

        ###NOTE consider adding a validation prompt.
    
    def plot_photo_fittings(self, trace_in_diary):
        photo_frames = get_frame_array(self.photo_470)

        fig, ax = plt.subplots(4)
        
        # plot 1: raw trace
        ax[0].plot(photo_frames, self.photo_470)
        ax[0].set_title('Raw 470')
        ax[0].set(ylabel=r'$\Delta$F/F')
        ax[0].set_xticklabels([])  # Remove x-tick labels
        
        # BIEXPONENTIAL FIT
        ax[1].set_ylim([self.photo_410.min(), self.photo_410.max()])
        ax[1].plot(photo_frames, self.photo_410)
        ax[1].plot(photo_frames, self.fit1, 'r')
        ax[1].set_title('410 with biexponential fit')
        ax[1].set(ylabel=r'$\Delta$F/F')
        ax[1].set_xticklabels([])  # Remove x-tick labels
        
        # 410 FIT TO 470
        ax[2].plot(photo_frames, self.photo_470)
        ax[2].plot(photo_frames, self.fit3, 'r')
        ax[2].set_title('410 fit over 470')
        ax[2].set(ylabel='F mean pixels')
        ax[2].set_xticklabels([])  # Remove x-tick labels
        
        # Bottom plot
        ax[3].plot(photo_frames, trace_in_diary)
        ax[3].set_xlabel('Frame #')  # x-label only on bottom plot
        ax[3].set_title('Smoothed trace')

        # Adjust spacing between subplots
        plt.subplots_adjust(hspace=0.5)  # Increase hspace for more vertical space between plots
        plt.savefig(f"{self.config.subject_path}\{self.config.trial_id}_fitting.png")
        
        if not self.config.silent_mode:
            plt.show(block=False)

        
    def identify_peak_positions(self, trace):
        # LIST OF PEAK POSITIONS
        x, _ = list(signal.find_peaks(trace, width=50, height=0.21))
        prom, l_base, r_base = list(signal.peak_prominences(trace, x, wlen=300))
        self.heights = trace[x]
        self.left_peak_bases = l_base
        return self.heights, self.left_peak_bases
    
    #@logger("log.txt")
    def get_clean_traces(self, photo_470, photo_410):
        
        self.fit_that_curve(photo_470, photo_410)

        trace_diary = self.normalize_and_zscore_trace(photo_470)
        trace_names = list(trace_diary.keys())

        #In case previous led correction was not successful
        self.plot_photo_fittings(trace_diary[trace_names[0]])
        if not self.config.silent_mode:
            while True:
                proper_led_assignment = input("Are 470 and 410 LEDs properly assigned? ('y'/ENTER or 'n'; input 's' to skip trial) \n NOTE: can silence this prompt by turning Silent mode on: ")
                plt.close()
                if proper_led_assignment.lower() in ['n', 'no']:
                    self.photo_410, self.photo_470 = photo_470, photo_410
                    self.fit_that_curve(photo_470, photo_410)
                    trace_diary = self.normalize_and_zscore_trace(photo_470)
                    self.plot_photo_fittings(trace_diary[trace_names[0]])
                    break

                if proper_led_assignment.lower() == 's':
                    self.config.report.append(self.trial_id)
                    print(f"Skipping {self.trial_id}")
                    return
                else:
                    break

            #In case isosbestic signal is not good.
            while True:
                good_isosbestic = input("Do you need to remove the isosbestic control signal? ('y'/'n'; press ENTER to skip): ")
                plt.close()
                if good_isosbestic.lower() in ['y', 'yes']: 
                    self.remove_isocesctic = True
                    trace_diary = self.normalize_and_zscore_trace(self.photo_470)
                    break
                else:
                    break       
        
        #get peak positions
        self.identify_peak_positions(trace_diary[trace_names[0]])

        return trace_diary 

class Restructure_Behavior_Data:
    def __init__(self, config, behav_raw, events_to_extract):
        """
        Description #NOTE
        
        Attributes:
        events_to_extract (lst of str): list of scored event information you wish to extract from photometry trace.
        event_dictionary (dict): Newly restructed event data that allows ethogram plotting
        """
        self.events_to_extract = events_to_extract
        self.behav_raw = behav_raw
        self.event_dictionary = {}
        self.zone = config.zone
        self.config = config
        self.add_zone = config.add_zone
        self.Groups = config.Groups
        self.trial_id = config.trial_id

    def create_dictionary(self):
        for event in self.events_to_extract:
            if event in self.config.events_to_score:
                bool_array = np.array(self.behav_raw[event], dtype = bool)
                self.event_dictionary.update({f"{event}": bool_array})

        return self.event_dictionary
    
    def define_zoned_behaviors(self, event_dictionary):
        if self.zone in self.events_to_extract:
            for event in self.events_to_extract:
                if event in [self.zone, self.config.start_parameter]:
                    continue
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
        else:
            raise ValueError(f"{self.config.zone} was not present in the 'Event extraction list'. Please be sure to include it when choosing events to analyze.")
                
        return self.event_dictonary

    #@logger("log.txt")
    def compile_behavior_groups(self):

        for new_grouping in self.Groups:
            compiled_array = np.zeros_like(self.behav_raw[self.config.start_parameter])
            for event in self.Groups[new_grouping]:
                if event not in self.behav_raw:
                    print(f"Event ({event}) to group into {new_grouping} is not present")
                    continue
                else:
                    compiled_array += self.behav_raw[event]
            clean_array = np.array(compiled_array, dtype = bool)
            self.event_dictionary.update({f'{new_grouping}' : clean_array})
        return self.event_dictionary
        
        
    def restructure(self):
        self.event_dictionary = self.create_dictionary()
        if self.Groups:
            self.event_dictionary = self.compile_behavior_groups()
        if self.add_zone:
            self.event_dictionary = self.define_zoned_behaviors()        
        return self.event_dictionary
            
class ExtractbyEvent:
    #@logger("log.txt")
    def __init__ (self, config, event_array, trace, time_markers):
        self.config = config
        self.event_array = event_array
        self.trace = trace
        self.offset = config.offset
        self.behav_fps = config.behav_fps
        self.photo_fps = config.photo_fps
        self.photo_convert = config.photo_fps / config.behav_fps
        self.start = time_markers[0]
        self.end = time_markers[1]
        self.photo_start = int(time_markers[0] * self.photo_convert)
        self.photo_end = int(time_markers[1] * self.photo_convert)
        self.min_duration = config.min_duration
        self.point_events = self.config.point_events
        self.left_peak_bases = []
        self.heights = []
        self.pre_s = config.pre_s
        self.post_s = config.post_s
        self.results = {}
        self.stop_at = config.perievent_limit
        self.Groups = config.Groups
        self.event_name = config.event_name
        self.UPPER_MARGIN = 0.1
        self.LOWER_MARGIN = 0.3
        self.font_size = 20
        self.subject_path = config.subject_path
        self.figure_path = make_folder('ROI_figures', self.subject_path)
        self.stage = config.stage
        self.photo_time = get_time_array(trace, self.photo_fps)
        self.trace_name = ''

    #@logger("log.txt")
    def extract_event_onsets(self, event_array, event_name):
        frame_onset = []

        for frame in range(self.start+1, self.end): #behavior params
            pre_frame = frame-1
            expected_min_duration = int(frame + floor((self.min_duration) * self.behav_fps))
            post_frame =frame+1

            try:
                if event_array[int(frame)] == 1 and event_array[pre_frame] == 0: #behavior
                    related_photo_frame = int(frame * self.photo_convert) #convert to photo
                    if event_name not in self.point_events:
                        if event_array[frame:expected_min_duration].all() == 1: #filters out potential score errors without impacting point events
                            frame_onset.append(related_photo_frame - (self.offset * self.photo_fps))
                        elif event_array[post_frame] == 1 and related_photo_frame in self.left_peak_bases: #prioritizes frame related to generated peak
                            frame_onset.append(related_photo_frame - (self.offset * self.photo_fps))
                    else:
                        frame_onset.append(related_photo_frame - (self.offset * self.photo_fps))
            except IndexError:
                print("error frame: ", frame, "value: ", event_name, "total frames: ", range(self.start+1, self.end), len(event_array)) 

        return frame_onset
        
    
    def clean_extracted_frames(self, list_of_frames):
        
        list1 = list_of_frames[1:]
        list2 = list_of_frames[:-1]

        for item1, item2 in zip(list1, list2):
            if int(item2) in range(item1, item1 + (self.pre_s * self.photo_fps)) or int(item2 - (self.pre_s - floor(self.pre_s/2)) * self.photo_fps) in range(item1, item1 + floor(self.post_s/1.5) * self.photo_fps): #1.5 is how much overlap is allowed.
                if item2 in list_of_frames:
                    list_of_frames.remove(item2)
        return list_of_frames
    
    def extract_photo_values(self, list_of_frames):
        """
        Splits the peri-event and baseline periods from the smooth trace.
        
        Parameters:
        smooth_trace (list or np.array): The smooth trace data.
        placement (list): The event placements.

        Returns:
        dict, dict, list, list: Peri-baseline, peri-event data, original placements, and adjusted placements.
        """
        peri_baseline = {} #peri-baseline has double the frames to match the event time in order to calculate area under the curve.
        peri_event = {}
        counter = 1

        # Iterate through each placement
        for i in list_of_frames:
            baseline_start = int(i - 1 - (self.post_s * self.photo_fps))  # Baseline start
            baseline_end = int(i - 1)  # Baseline end
            event_end = int(i + (self.post_s * self.photo_fps))      # Event end
            
            # Only process valid ranges within the smooth trace boundaries
            if baseline_start >= 0 and event_end + 5 <= len(self.trace):
                # Retrieve baseline and event values using list comprehension
                baseline = self.trace[baseline_start:baseline_end]
                event = self.trace[i:event_end]

                # Store in the respective dictionaries
                peri_baseline[counter] = baseline
                peri_event[counter] = event
                counter += 1

        return peri_baseline, peri_event, list_of_frames       
     
    #@logger("log.txt")
    def clean_restructure_photo_values(self, peri_baseline, peri_event, frame_onset, event):  #good_values
        """
        Compute the difference between baseline and event means, filtering by 'stop_at'.
        The results are returned in a dictionary to avoid verbosity.

        Parameters:
        peri_baseline (dict): Baseline data.
        peri_event (dict): Event data.
        placement (list): Placement data.
        i (int): Group identifier.

        Returns:
        dict: Contains event means, baseline means, baseline matrix, event matrix, placements, and diff.
        """
        # Compute means for baseline (cropped) and event
        event_means = [mean(peri_event[key]) for key in peri_event]
        baseline_means = [mean(peri_baseline[key][-int(self.pre_s*self.photo_fps):]) for key in peri_baseline]
        
        # Calculate percentage differences
        diff = [( (b - a) / a ) * 100 for a, b in zip(baseline_means, event_means)]

        result = {} 

        if self.stop_at > 0:
            try:
                # Create and sort DataFrame
                workme = pd.DataFrame({
                    'event means': event_means,
                    'baseline means': baseline_means,
                    'placement_s': (frame_onset * self.photo_convert),
                    'baseline_matrix': list(peri_baseline),
                    'event_matrix': list(peri_event),
                    'diff': diff
                }).sort_values(by='diff', ascending=False)

                # Select top results based on `stop_at`
                top_results = int(self.stop_at * 2) if event in self.Groups.keys() else int(self.stop_at)
                df = workme.head(top_results)

                # Update dictionary with sorted values
                result.update({
                    "event_means": df['event'].tolist(),
                    "baseline_means": df['baseline'].tolist(),
                    "placement_s": df['placement'].tolist(),
                    "peri_baseline_matrix": np.array([peri_baseline[key] for key in df['base2']], dtype=float),
                    "peri_event_matrix": np.array([peri_event[key] for key in df['event2']], dtype=float)
                })
            
            except ValueError:
                # If there's an issue, return empty results
                print("Data Error: Returning empty values")
                result.update({
                    "event_means": [],
                    "baseline_means": [],
                    "placement_s": [],
                    "peri_baseline_matrix": np.array([]),
                    "peri_event_matrix": np.array([])
                })

        else:
            # Default case without filtering
            result.update({
                'event_means': event_means,
                'baseline_means': baseline_means,
                'placement_s': [i * self.photo_convert for i in frame_onset],
                "peri_baseline_matrix": np.array([peri_baseline[key] for key in peri_baseline], dtype=float),
                "peri_event_matrix": np.array([peri_event[key] for key in peri_event], dtype=float)
            })
        
        return result
    
    def rep_plot_overlap(self):

        if len(self.event_array) != 0:
            restructured_frame_onset = []
            frame_end = []

            for frame in range(1, len(self.event_array)):
                if self.event_array[frame - 1] == False and self.event_array[frame] == True:
                    restructured_frame_onset.append(floor(frame * self.photo_convert))
                if self.event_array[frame] == True and self.event_array[frame + 1] == False:
                    frame_end.append(floor(frame * self.photo_convert))
                if self.event_array[-2] == True and self.event_array[-1] == True:
                    frame_end.append(floor(frame * self.photo_convert))

            plt.figure(figsize = (10, 5))

            plt.plot(self.photo_time, self.trace, color = 'k')
            plt.ylim([self.trace.min() - (self.UPPER_MARGIN * self.trace.max()), self.trace.max() + (self.LOWER_MARGIN * self.trace.max())])
            plt.ylabel(r'zF', fontsize = self.font_size)
            plt.xlabel('Time (s)', fontsize = self.font_size)

            plt.xticks(fontsize = self.font_size); plt.yticks(fontsize = self.font_size)

            plt.tight_layout()

            plt.axvspan(self.photo_time[self.photo_start], self.photo_time[self.photo_end], color = 'palegoldenrod', alpha = 0.2, lw = 0)
            for roi in restructured_frame_onset:
                if self.photo_start <= roi <=self.photo_end:
                    postevent = int(roi + self.config.post_s * self.config.photo_fps)
                    preevent = int(roi - self.config.pre_s * self.config.photo_fps)
                    plt.axvspan(self.photo_time[roi], self.photo_time[postevent], color = 'coral', alpha = 0.5, ymax = 0.6)
                    plt.axvspan(self.photo_time[roi], self.photo_time[preevent], color = 'aquamarine', alpha = 0.5, ymax = 0.6)
            for s, e in zip(restructured_frame_onset, frame_end):
                try:
                    plt.axvspan(self.photo_time[s], self.photo_time[e], color = 'midnightblue', ymax = 0.1)
                except IndexError:
                    continue
            plt.title(f"{self.event_name}")
            plt.savefig(f'{self.figure_path}\{self.event_name}_{self.stage}_{self.config.trial_id}')
            plt.close('all')

    # def photo_ethogram_plot(self, event_dictionary, color_dictionary):
        

    def get_values(self):
        frame_onset = self.extract_event_onsets(self.event_array, self.event_name)
        clean_frame_onset = self.clean_extracted_frames(frame_onset)
        peri_baseline, peri_event, frame_onset = self.extract_photo_values(clean_frame_onset)
        return self.clean_restructure_photo_values(peri_baseline, peri_event, frame_onset, self.event_name)

def save_variables(result, config, behavior_path):

    for key in result:
        if len(result[key]) != 0:
            np.savetxt(f'{behavior_path}/{key}_{config.trial_id}.csv', result[key], delimiter = ',', fmt = '%s')

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
    return np.array([i for i in range(len(array1))])

def get_time_array(array2, fps):
    return np.array([i/fps for i in range(len(array2))])

#@logger("log.txt")
def make_folder(new_folder, parent_directory):
    """
    Creates new folders and or sets the directory. 

    Args:
        parent_directory (str): The parent directory for the new folder.
        new_folder (str): The name of the new fodler to be created.

    Returns:
        full_path (str): The new directory where the folder was created
    
    Raises: 
        FileNotFoundError: If the specified parent directory does not exist.
        PermissionError: If the directory is not writable or the folder cannot be created.

    """

    mode = 0o666

    full_path = os.path.join(parent_directory, new_folder)

    #check if parent directory exists
    if not os.path.exists(parent_directory):
        raise FileNotFoundError(f"Parent directory '{parent_directory}' does not exist.")
    
    #checks if user has permission to write to that directory
    if not os.access(parent_directory, os.W_OK):
        raise PermissionError(f"Write permission denied for directory '{parent_directory}.")
    
    #Creates the folder if it doesnt exists
    if not os.path.exists(full_path):
        try:
            os.mkdir(full_path, mode)
        except OSError:
            raise PermissionError(f"Failed to create directory {full_path}. Check permissions: {OSError}")


    return full_path

def exp2(x, a, b, c, d):
    return a*exp(b*x) + c*exp(d*x)

#@logger("log.txt")
def main(configuration_file):
    config = AnalysisParameters(load_config(configuration_file))

    ##NOTE Prompt for GUI to select behaviors to extract

    #Begin iterating through analysis folders
    analysis_folders = os.listdir(config.project_home)
    event_save_path = make_folder("Behaviors", config.project_home)
    for subject_trial_id in analysis_folders:
        #Exclude processed value or video folders that may be in the project home directory
        config.trial_id = subject_trial_id
        exclude_folder = ['Behaviors', 'Videos', 'Summary', 'Archive']
        if not os.path.isdir(os.path.join(config.project_home, subject_trial_id)) or subject_trial_id in exclude_folder:
            continue
        else:
            #Identify and load raw files
                      
            processor = Clean_and_ProcessFiles(config); 
            photo_470, photo_410, events_to_extract, behav_raw_df = processor.process_files()
            
            config.subject_path = processor.subject_trial_path
            
            plotter = SmoothPlotnValidate_Photo(config, photo_470, photo_410)
            trace_diary = plotter.get_clean_traces(photo_470, photo_410)
            if trace_diary is None:
                continue

            redo_behavior = Restructure_Behavior_Data(config, behav_raw_df, events_to_extract)
            event_dictionary = redo_behavior.restructure()

            if config.controls:
                stages = ['before', 'control', 'during']
            else:
                stages = ['before', 'during']
            for trace in trace_diary:
                trace_path = make_folder(trace, event_save_path)
                for event_name in event_dictionary:
                    config.event_name = event_name
                    for stage in stages:
                        stage_path = make_folder(stage, trace_path)
                        # print(stage, event_name)
                        config.stage = stage
                        try:
                            if stage == 'before':
                                time_markers = [0, processor.start_times[0]]
                            elif stage == 'control':
                                time_markers = [processor.start_times[0], processor.end_times[0]]
                            elif stage == 'during' and config.controls:
                                time_markers = [processor.start_times[1], processor.end_times[1]]
                            elif stage == 'during' and not config.controls:
                                time_markers = [processor.start_times[0], processor.end_times[0]]
                        except IndexError:
                            print("Check 'use_control_trial' parameter in configuration file.")

                        event_path = make_folder(event_name, stage_path)
                        extraction = ExtractbyEvent(config, event_dictionary[event_name], trace_diary[trace],  time_markers)
                        results = extraction.get_values()
                        save_variables(results, config, event_path)
                # plot_photometry_ethogram(config, event_dictionary[event_name], trace, time_markers)
                        if trace == list(trace_diary.keys())[0]:
                            extraction.rep_plot_overlap()


if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(configuration_file)
