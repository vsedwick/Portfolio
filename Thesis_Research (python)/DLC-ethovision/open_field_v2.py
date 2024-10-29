###last updated 10/29/2024

"""
This script is meant to post-process open field videos that were labeled using deeplabcut.
The purpose is to calculate distance traveled, time spent in specific areas, 
and plot the tracklets within the arena. Additionally, this script will clean and smooth
aberrant tracklets and establish arena zones. 

This script was written to accompany the OpenField DLC profile trained by Victoria Sedwick in the Autry Lab.

If you have any questions or concerns, please contact Victoria Sedwick at sedwick.victoria@gmail.com.
"""

config_file = "config.yaml"  ##NOTE EDIT

import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon
import os
import matplotlib.pyplot as plt
import sys
import math
from datetime import datetime
import yaml
from ruamel.yaml import YAML


def read_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)

    except FileNotFoundError:
        try:
            config_path = os.path.join(os.getcwd(), config_file)
            with open(config_path, 'r') as file:
                config = yaml.safe_load(file)        
        except FileNotFoundError: 
            print(f"Error: Configuration file {config_file} not found.")
            sys.exit(1) 
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
        sys.exit(1)


    return config

def markers( x):
    """
    Identifies the bodypart markers that are in each video DataFrame
    """
    bodyparts = [*set(list({i[0] for i in x.columns.values}))]
    return bodyparts

def load_data(file):
    dataframe = pd.read_hdf(file)
    dataframe.head(); scorer = dataframe.columns.get_level_values(0)[0]
    dataframe = dataframe[scorer]
    labels = markers(dataframe)

    return dataframe, labels

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

def update_config(configuration_params, config_file):
    yaml = YAML()
    yaml.preserve_quotes = True  # Optional: Preserves quotes around values

    # Load the YAML file with comments
    with open(config_file, 'r') as file:
        config_data = yaml.load(file)
    
    # Update root-level configurations
    config_data['scoretype'] = configuration_params.scoretype  # Corrected typo from 'scoretpye'

    # Update Video Parameters
    video_params = config_data.get('Video Parameters', {})
    video_params['real_sizecm'] = configuration_params.real_sizecm
    video_params['length_s'] = configuration_params.length_s
    config_data['Video Parameters'] = video_params
    #NOTE add mcf validation and clear the variable

    # Update Arena Parameters
    arena_params = config_data.get('Arena_parameters', {})
    arena_params['scale_for_center_zone'] = configuration_params.center_scale
    arena_params['plikelihood'] = configuration_params.likelihood
    arena_params['bodyparts_to_analyze'] = configuration_params.trackBP
    arena_params['bodyparts_to_determine_trial_start'] = configuration_params.startBP
    config_data['Arena_parameters'] = arena_params

    # Write the updated data back, preserving comments
    with open(config_file, 'w') as file:
        yaml.dump(config_data, file)

class VideoAnalysis:

    def __init__(self, config):

        self.config_params = self.ConfigHandler(config)

        self.trackBP = self.config_params.trackBP
        self.config = config
        self.Videos = []; self.Distance=[]; self.Velocity=[]; self.Time_in_Center=[]
        self.Time_in_Border=[]; self.Time_in_Arena=[]; self.BorderEntry=[]; 
        self.Time_immobile = []; self.Time_in_motion = []; self.CenterEntry=[]
        self.rootdir = config["rootdir"]

        self.file_name = None

        self.image_save = make_folder("Tracklets_plots", self.rootdir)

    def get_timestamp(self):
        """
        Grabs current timestamp
        """
        dt = str(datetime.now()).replace(' ', '_').replace(':', '-').split('.')[0]
        return dt
    
    def video_name(self, file):
        """
        Extracts video name for final dataframe output
        """
        if '/' in file:
            file = file.split('/')[-1]
        if '\\\\' in file:
            file = file.split('\\\\')[-1]
        if '\\' in file:
            file = file.split('\\')[-1]
        return file.replace('DLC', '#').split('#')[0]


    def analyze(self, file, show_plots = True):
        
        self.file_name = self.video_name(file)

        self.Videos.append(self.file_name )

        #load file and configuration paramters
        self.config_params.dataframe, self.config_params.labels = load_data(file)

        #initialize Arena param

        arena_params = self.Arena(self.config_params)
        arena_coords, center_coords, wall_coords = arena_params.make_arena()

        plt.plot(arena_coords[0], arena_coords[1], 'k')
        plt.plot(center_coords[0],center_coords[1], 'k')
        plt.plot(wall_coords[0], wall_coords[1], 'k')

        if self.config_params.length_s != 0:
            start, end = arena_params.crop_video()
            self.config_params.dataframe = self.config_params.dataframe[start:end]
        
        for bp in self.config_params.trackBP:
            x, y = self.config_params.dataframe[bp]['x'].values, self.config_params.dataframe[bp]['y'].values

            get_metrics = self.Metrics(self.config_params, arena_params, [x, y], bp)

            #clean and interpolate tracklets
            x, y = get_metrics.cleanup()
            x, y = get_metrics.fill_in()

            #calculate distance and velocity
            dist_vel, motion, time_zone, entries = get_metrics.run_calculation()

            self.Distance.append(dist_vel[0])
            self.Velocity.append(dist_vel[1])
                    
            self.Time_immobile.append(motion[1])
            self.Time_in_motion.append(motion[0])
                    

            self.Time_in_Center.append(time_zone[1])
            self.Time_in_Border.append(time_zone[2])
            self.Time_in_Arena.append(time_zone[0])

            #Frequency in zones#

            self.BorderEntry.append(entries[0])
            self.CenterEntry.append(entries[1])

            plt.plot(x,y, label=f'{bp}')
            
        plt.legend()
        if show_plots:
            plt.show(block = False)
        if self.config_params.save_images:
            image_path = os.path.join(self.image_save, f"{self.file_name}.png")

            plt.savefig(image_path)


    def save_df(self):
        
        summary_dict = {
                        'Total Distance Traveled (cm)': self.Distance, 
                        'Velocity (cm/s)': self.Velocity, 
                        'Center Entries': self.CenterEntry,
                        'Time in Center': self.Time_in_Center,
                        'Border Entry': self.BorderEntry,
                        'Time in Border': self.Time_in_Border,
                        'Time in Arena': self.Time_in_Arena, 
                        "Time Immobile": self.Time_immobile,
                        "Time Moving": self.Time_in_motion
                        }
        # print(summary_dict)
        # exit()
        
        if self.config_params.scoretype.lower() == 'single':
            video_summary=pd.DataFrame(summary_dict, index = self.config_params.trackBP)
            summary=os.path.join(self.rootdir, f'{self.file_name}-{self.get_timestamp()}.xlsx')
            video_summary.to_excel(summary, index=True, header=True)

            print(f"Analysis for {self.file_name} is saved")
            sys.exit()

        else:
            headersss=[]
            headersss.append(self.Videos)
            headersss.append(self.config_params.trackBP)
            ind = pd.MultiIndex.from_product(headersss, names=["video", "label"])
            video_summary=pd.DataFrame(summary_dict, index = ind)

            summary=os.path.join(self.rootdir, f'{self.config_params.project_name}-{self.get_timestamp()}.xlsx')
            video_summary.to_excel(summary, index=True, header=True)

            print(f"{self.config_params.project_name} is saved.")

    class ConfigHandler:
        def __init__(self, config):
            
            self.calibration_mode = config["calibration_mode"]

            self.save_images = config["Arena_parameters"]["save_tracklets"]
            
            self.examp_file = config["example_behavior_file_h5"]

            self.file = None

            self.dataframe, self.labels = load_data(self.examp_file)

            self.project_name = config["project_name"]

            self.real_sizecm = self._val_real_size(config["Video Parameters"]["real_sizecm"])
            self.fps = self._val_fps(config["Video Parameters"]["fps"])
            
            self.mcf = config["Arena_parameters"]["movement_correction_factor"]
            self.center_scale = self.val_center(config["Arena_parameters"]["scale_for_center_zone"])
            
            self.trackBP = []
            self.startBP = []
            self.trackBP = self._validate_bp(config.get("Arena_parameters", {}).get("bodyparts_to_analyze", []), "Analysis")
            
            print("trackBP", self.trackBP)
            self.length_s = self._val_crop_length(config["Video Parameters"]["length_s"])
            self.startBP = self._validate_bp(config.get("Arena_parameters", {}).get("bodyparts_to_determine_trial_start", []),  "Starter")
            self.scoretype = self._val_scoretype(config["scoretype"])
            self.likelihood = self._val_likelihood(config["Arena_parameters"]["plikelihood"])


        def _validate_bp(self, lst, cat):
            #NOTE does not handle mispellings properly
            
            if isinstance(lst, list):
                for i in range(len(lst)):
                    if lst[i] not in self.labels:
                        if cat == "Analysis":
                            print(f"BP {lst[i]} in variable 'bodyparts_to_analyze' is not present in your dataset")
                        if cat == "Starter":
                            print(f"BP {lst[i]} in variable 'bodyparts_to_determine_trial_start' is not present in your dataset")
                        lst[i] == None
                new_list = [i for i in lst if i is not None]  
                return new_list
                
            elif lst is None:
                if cat == "Analysis":
                    self.trackBP = self.track_labels()
                    print("deed is done")
                    return self.trackBP

                if cat == "Starter":
                    self.startBP = self.start_bp()
                    return self.startBP

            else:
                return []
            
        def _val_scoretype(self, scoretype):
            # Determine score type
            options = ['single', 'batch', 'multiple']
            while not scoretype or scoretype not in options:
                new_scoretype = input("Are you analyzing a single (type 'single') or multiple (type 'batch' or 'multiple') videos?: ").lower()
                if new_scoretype in options:
                    return new_scoretype
                else:
                    print(f"Invalid input. Please enter {options}.")
                    continue
            return scoretype

        def _val_real_size(self, real_sizecm):
            if not real_sizecm or real_sizecm < 30:
                while True:
                    new_real_sizecm = int(input("What is the size of the arena (hypotenuse in centimeters e.g. 65): "))
                    if isinstance(new_real_sizecm, int) and new_real_sizecm > 30:
                        return new_real_sizecm
            else:
                return real_sizecm
        def _val_fps(self, fps):
            if not fps:
                while True:
                    new_fps = int(input("What is the frame rate of the videos to be analyzed?: "))
                    if isinstance(new_fps, int) and self.real_sizecm > 5:
                        return new_fps
            else:
                return fps
        def val_center(self, center_scale):
            # Set center scale
            if not center_scale or 0 > center_scale > 1:
                while True:
                    try:
                        scale = float(input("Set scale for center zone between 0 and 1 (0.28 is most similar to Ethovision): "))
                        if 0 < scale < 1:
                            self.center_scale = scale
                            return self.center_scale
                        else:
                            print('Scale is outside of range. Please enter a value between 0 and 1.')
                    except ValueError:
                        print("Invalid input. Please enter a numeric value between 0 and 1.")
            else:
                return center_scale
            
        def _val_crop_length(self, length_s):
        # Crop videos if needed
            if length_s is None: 
                while True:
                    crop = input("Do you want to crop your videos? (yes or no); Start time will be determined by a certain likelihood of the mouse nose: ").lower()
                    if crop in ['yes', 'y']:
                        while True:
                            try:
                                new_length_s = int(input("How much time do you wish to analyze in seconds?: "))
                                if new_length_s > 0:
                                    return new_length_s
                                else:
                                    print("Value needs to be greater than 0.")
                            except ValueError:
                                print("Invalid input. Please enter a positive integer.")
                    elif self.crop == 'no':
                        return 0
                    else:
                        print("Invalid input. Please enter 'yes' or 'no'.")
            else:
                return length_s
            
        def _val_likelihood(self, likelihood):        # Determine likelihood of accuracy
            if not likelihood:
                while True:
                    try:
                        new_likelihood = float(input('What is the likelihood cutoff e.g. >0.8. (Input number <1): '))
                        if 0 < new_likelihood < 1:
                            self.likelihood = new_likelihood
                            return self.likelihood
                        else:
                            print("Likelihood must be a number between 0 and 1.")
                    except ValueError:
                        print("Invalid input. Please enter a number between 0 and 1.")
            return likelihood
        
        def start_bp(self):
            print(self.labels)
            k = 0
            
            while k < 3:
                sBP = str(input("\n What bodyparts (3) should be in the arena at the start of tracking? If only one, enter 3 times: "))
                if str(sBP) in self.labels:
                    k += 1
                    self.startBP.append(str(sBP))
                else:
                    print("That bodypart is not in this file")
                    continue
            
            print(self.startBP)
            return self.startBP

        def track_labels(self):
            #NOTE should add check for 'in list'
            # Adding labels
            while True:
                print(self.labels)
                label = input("What body marker would you like to track?: ")
                if label in self.labels:
                    self.trackBP.append(label) 
                    break

            while True:
                more = input("Add another? Input label or type 'no' or ENTER to skip: ")
                if more.lower() in ['no', 'n', '']:
                    break
                elif more in self.labels:
                    self.trackBP.append(more)
                    print(self.trackBP)
                else:
                    print(f"Label '{more}' not found in the available labels. Try again.")

            # Removing labels
            while True:
                bye = input("Remove label(s)? Input label or type 'no' or ENTER to skip: ")
                if bye.lower() in ['no', 'n', '']:
                    break
                elif bye in self.trackBP:
                    self.trackBP.remove(bye)
                    print(self.trackBP)
                else:
                    print(f"Label '{bye}' is not in the current list. Try again.")

            return self.trackBP
        
        def _val_h5(self, file):
            if not file.endswith('.h5'):
                raise ValueError("Please re-load 'example_behavior_file_h5' file with a '.h5' extension")
            else:
                return file

        def clear_variables(self):
            self.length_s = self.real_sizecm = self.likelihood = self.startBP = self.center_scale = self.trackBP = None


    
    class Arena:
        def __init__(self, config_params):
            self.config_params = config_params
            self.dataframe = config_params.dataframe
            self.arena = None
            self.center_zone = None
            self.border = None
            self.border1 = None
            self.wall = None
            self.Frame = None
            self.startBP = config_params.startBP
            self.trackBP = config_params.trackBP
            self.center_scale = config_params.center_scale
             
        def make_arena(self):
            """
            Sets the dimensions of the open field arean based on the arena labels in the video.
            The center zone is scaled using the "center_scale" parameter.
            """
            #Index the x coordinates of the label and average them all for a static point
            TLx = self.dataframe['TL']['x'].mean()
            TLy = self.dataframe['TL']['y'].mean()
            TRx = self.dataframe['TR']['x'].mean()
            TRy = self.dataframe['TR']['y'].mean()
            BLx = self.dataframe['BL']['x'].mean()
            BLy = self.dataframe['BL']['y'].mean()
            BRx = self.dataframe['BR']['x'].mean()
            BRy = self.dataframe['BR']['y'].mean()

            #outter arena
            self.arena = np.array([[TLx, TLy], [TRx, TRy], [BRx, BRy], [BLx, BLy]])

            #center zone
            difference = self.scale_arena()
            grow = self.grow_arena()
            self.arena = Polygon(self.arena)
            self.center_zone = Polygon(self.arena.buffer(-difference))
            self.arena1 = Polygon(self.arena.buffer(+grow, cap_style = 3))
            
            #establish borders
            center_zonex,center_zoney=self.center_zone.exterior.xy
            arenax,arenay=self.arena.exterior.xy

            bordery=center_zoney+arenay
            borderx=center_zonex+ arenax

            border=[[borderx[i], bordery[i]] for i in range(len(borderx))] 
            self.border1=Polygon(border)

            arenax,arenay=self.arena.exterior.xy

            bordery=center_zoney +arenay
            borderx=center_zonex + arenax

            border=[[borderx[i], bordery[i]] for i in range(len(borderx))] 
            self.border=Polygon(border)

            #establish walls
            arena1x,arena1y = self.arena1.exterior.xy
            borderx, bordery = self.border.exterior.xy

            wally = arenay + arena1y
            wallx = arenax + arena1x

            wall=[[wallx[i], wally[i]] for i in range(len(wallx))] 
            self.wall = Polygon(wall)

            return [arenax, arenay], [center_zonex, center_zoney], [wallx, wally]

    
        def scale_arena(self):
            """
            Creates the center zone
            """
            xs = [i[0] for i in self.arena]
            ys = [i[1] for i in self.arena]
            x_center = 0.5 * (min(xs) + max(xs))
            y_center = 0.5 * (min(ys) + max(ys))
            center = Point(x_center, y_center)
            min_corner = Point(min(xs), min(ys))
            return center.distance(min_corner) * self.center_scale

        def grow_arena(self):
            """Creates the "walls" to account for tracklets outside of the arena floor such as when the subject is rearing
            """

            xs = [i[0] for i in self.arena]
            ys = [i[1] for i in self.arena]
            x_center = 0.5 * (min(xs) + max(xs))
            y_center = 0.5 * (min(ys) + max(ys))
            center = Point(x_center, y_center)
            min_corner = Point(min(xs), min(ys))
            return center.distance(min_corner) * 0.15
        


        def crop_video(self):
            """
            crops video based on user input
            
            """

            k=0
            crop=[]
            standard1=self.dataframe[self.startBP[0]]['likelihood'].values
            standard2=self.dataframe[self.startBP[1]]['likelihood'].values
            standard3=self.dataframe[self.startBP[2]]['likelihood'].values
            standard=tuple(zip(standard1, standard2, standard3))
            for (i, j,) in zip(standard, range(len(standard))):
                a,b, c=i
                if a>=0.9 and b>=0.9 and c>=0.9:
                    k+=1
                else:
                    k=0
                if k==3:
                    crop.append(j-2)
                    break
                else:
                    continue
            start=crop[0]
            crop.append((start+((self.config_params.length_s)*self.config_params.fps)))
            end=crop[1]
            return start, end;

        def plot_arena_template(self):
            self.make_arena()
            arenax, arenay = self.arena.exterior.xy
            centerx, centery = self.center_zone.exterior.xy
            wallx, wally = self.wall.exterior.xy
            borderx, bordery = self.border.exterior.xy
            
            plt.figure("Arena")
            plt.plot(arenax, arenay, 'k')
            plt.fill(arenax, arenay, 'red', alpha=0.4, label='Arena')
            plt.plot(centerx, centery, 'k')
            plt.fill(centerx, centery, 'y', alpha=0.3, label='Center')
            plt.plot(wallx, wally, 'k')
            plt.fill(wallx, wally, 'lime', alpha=0.3, label='Wall')
            plt.plot(borderx, bordery, 'k')
            plt.fill(borderx, bordery, 'b', alpha=0.3, label='Borders')
            plt.axis('equal')
            plt.legend(loc='lower left')
            plt.show(block=False)
            
        def validate_arena(self):
            
            # Create and validate arena
            while True:
                self.plot_arena_template()

                good = input("Is this arena okay? yes(or ENTER key) or no: ").lower()
                plt.close('all')
                if good in ['yes', 'y', '']:
                    break
                elif good in ['no', 'n']:
                    print("Try adjusting the re-adjusting the center scale. If the problem persist, review your 'arena markings' (e.g. 'TL', 'TR', etc.)")
                    self.config_params.center_scale = self.config_params.val_center(False)
                    continue

    class Metrics:
        def __init__(self, config_params, arena_params, coords, bp):
            self.arena1 = arena_params.arena1
            self.arena = arena_params.arena
            self.center_zone = arena_params.center_zone
            self.border = arena_params.border

            self.dataframe = config_params.dataframe
            
            self.x = coords[0]
            self.y = coords[1]

            self.likelihood = config_params.likelihood
            self.real_sizecm = config_params.real_sizecm
            self.fps = config_params.fps
            self.mcf = config_params.mcf
            self.bp = bp

        def cleanup(self):
            coords = list(zip(self.x, self.y))
            cleaned_coords = [(a, b) if Point(a, b).within(self.arena1) else (np.nan, np.nan) for a, b in coords]
            x_clean, y_clean = zip(*cleaned_coords)
            x_series = pd.Series(x_clean).interpolate(method='linear')
            y_series = pd.Series(y_clean).interpolate(method='linear')
            return x_series, y_series

        def fill_in(self):
            plikelihood = self.dataframe[self.bp]['likelihood'].values
            coords = list(zip(self.x, self.y))
            filled_coords = [(a, b) if lik >= self.likelihood else (np.nan, np.nan) for (a, b), lik in zip(coords, plikelihood)]
            x_filled, y_filled = zip(*filled_coords)
            x_series = pd.Series(x_filled).interpolate(method='linear')
            y_series = pd.Series(y_filled).interpolate(method='linear')
            return x_series, y_series

        def calculate_distance(self):
            arenax, arenay = self.arena.exterior.xy
            dist_arena = float(math.sqrt((arenax[2] - arenax[0]) ** 2 + (arenay[2] - arenay[0]) ** 2))
            real_factor = self.real_sizecm / dist_arena
            coords = list(zip(self.x[:-1], self.y[:-1], self.x[1:], self.y[1:]))
            distances = [float(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)) * real_factor * self.mcf for x1, y1, x2, y2 in coords]
            return sum(distances)

        def calculate_velocity(self, distance):
            return distance / (len(self.x)/self.fps)

        def calculate_movement(self):
            coords = list(zip(self.x[:-1], self.y[:-1], self.x[1:], self.y[1:]))
            move_count = sum(1 for x1, y1, x2, y2 in coords if float(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)) > 0.8)
            no_move_count = len(coords) - move_count
            return move_count / self.fps, no_move_count / self.fps

        def time_in_zones(self):
            coords = list(zip(self.x, self.y))
            center_time = sum(1 for a, b in coords if Point(a, b).within(self.center_zone)) / self.fps
            arena_time = sum(1 for a, b in coords if Point(a, b).within(self.arena1)) / self.fps
            border_time = sum(1 for a, b in coords if Point(a, b).within(self.border)) / self.fps
            return arena_time, center_time, border_time

        def frequency_crossing(self):
            coords = list(zip(self.x[:-1], self.y[:-1], self.x[1:], self.y[1:]))
            border_entry = sum(1 for x1, y1, x2, y2 in coords if Point(x1, y1).within(self.center_zone) and Point(x2, y2).within(self.border))
            center_entry = sum(1 for x1, y1, x2, y2 in coords if Point(x1, y1).within(self.border) and Point(x2, y2).within(self.center_zone))
            return border_entry, center_entry

        def run_calculation(self):
            dist = self.calculate_distance()
            vel = self.calculate_velocity(dist)
            dist_vel = [dist, vel]

            motion = self.calculate_movement()

            time_zone =  self.time_in_zones()

            entries = self.frequency_crossing()



            return dist_vel, motion, time_zone, entries

def main(config_file):

    #load directory information
    config = read_config(config_file)    

    #confirm or adjust the arena parameters
    calibration_mode = config.get("calibration_mode", True)

    if calibration_mode:
        
        while True:
            config = read_config(config_file) 
            #will update and adjust missing variables
            calibrator = VideoAnalysis(config)

            arena_maker = calibrator.Arena(calibrator.config_params)
            #confirm Arena is right
            arena_maker.validate_arena()

            #finish analysis for example file
            calibrator.analyze(calibrator.config_params.examp_file, show_plots = True)
            #validate calibration
            confirmation = input("Does everything look okay? \n Input 'yes' or 'y' to continue with analysis or 'no' or 'restart' to change parameters: ")
            if confirmation.lower() in ['yes', 'y', '']:
                #need to update config file.
                plt.close()
                update_config(calibrator.config_params, config_file)
                #save info
                if calibrator.config_params.scoretype.lower() == 'single':
                    calibrator.save_df()
                    exit('"Single" analysis complete.')
                
                break
            
            else:
                clear = input("Would you like to clear your configuration file and restart the setup?: ")
                if clear.lower() in ['yes', 'y']:
                    calibrator.config_params.clear_variables()
                    update_config(calibrator.config_params, config_file)
                    continue
                else:
                    exit("Exiting program")

    print("Beginning batch analysis, please wait... ")

    final_config = read_config(config_file)
    batch_analyzer = VideoAnalysis(final_config)
    all_files = os.listdir(batch_analyzer.rootdir)
    file_list = [i for i in all_files if i.endswith('.h5')]

    for file in file_list:
        file_path = os.path.join(batch_analyzer.rootdir, file)
        batch_analyzer.analyze(file_path, show_plots = False)
        plt.close()
    
    batch_analyzer.save_df()

    print(f"Analysis complete. Tracklets can be found in {batch_analyzer.image_save}")



if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(config_file)
