
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

    ##Prompt for GUI to select behaviors to extract

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
        self.behavior_time = []
        self.behav_frames = []
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

    def get_time_and_frames(self):

        self.behav_frame = [i for i in range(len(self.behav_raw))]
        self.behav_time = [i/self.config.behav_fps for i in self.behav_frame]
    
        self.photo_frame = [i for i in range(len(self.photo_470))]
        self.photo_ime = [i/self.config.photo_fps for i in range(len(self.photo_470))]   
    
    def validate_files(self):
        if self.config.start_parameter not in self.behaviors:
            self.config.report.append(self.subject_trial_id)
            return False
        else:
            return True
    
    def set_times_equal(self): #might have 1 sec difference
        if self.behavior_time.max()+1>self.photo_time.max():
            cutit = floor(self.config.offset)

            j = self.behav_raw[self.behaviors[0]]
            self.behav_frames = np.array([i for i in range(len(j[cutit:]))])
            self.behavior_time = np.array([i/self.config.behav_fps for i in self.behav_frames])

            self.behav_raw = self.behav_raw[cutit:]
        
        return self.behavior_time, self.behav_raw

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

            #Generate time and frame arrays? Maybe.
        else:
            raise ValueError("Check the 'full_trace' paramter in the configuration file. Make sure it says 'yes'.") 
        
    def process_files(self):
            self.open_files()
            self.split_behavior_values()

            if not self.validate_files():
                return False
            else:
                self.split_photometry_values()

                self.get_time_and_frames()
                #Correct frame rate
                print(f"Correcting frames per second. Current fps: {float(len(self.photo_470)/self.behavior_time.max())}")
                self.correct_fps()

                self.photo_time = [i/self.config.photo_fps for i in len(self.photo_470)]
                print(f"Corrected frames per second: {len(self.photo_470)/max(self.photo_time)}")
                
                #prepare to crop
                print(f"Time difference between video and trace before cropping: {self.behav_time.max()}, {self.photo_time.max()}")
                if self.config.full_trace == 'yes':
                    self.crop_file_ends()
                else:
                    self.crop_around_start_param()
                
                print(f"Time difference between video and trace after cropping: {self.behav_time.max()}, {self.photo_time.max()}")
                
                self.get_time_and_frames()

                return self.photo_410, self.photo_410, self.behaviors, self.behav_raw
    
    def correct_fps(self):
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

class SmoothNPlot:
    def __init__(self, config, photo_470, photo_410):
        self.photo_470 = self.photo_470
        self.photo_410
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



if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(configuration_file)


        

