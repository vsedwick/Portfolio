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

def main(config_file):

    config, rootdir, example_behav_file, file_list = read_config(config_file)    

    #confirm or adjust the arena parameters
    if example_behav_file.endswith('.h5'):
        calibrator = VideoAnalysis(config)
        while True:
            calibrator.set_up_analysis(example_behav_file)
            confirmation = input("Does everything look okay? \n Input 'yes' or 'y' to continue with analysis or 'no' or 'restart' to change parameters")
            if confirmation.lower() in ['yes', 'y', '']:
                break
            else:
                clear = input("Would you like to clear your configuration file and restart the setup?: ")
                if clear.lower() in ['yes', 'y']:
                    calibrator.clear_variables()
                    continue
                else:
                    exit("Exiting program")
    else:
        exit("Please re-load 'example_behavior_file_h5' file with a '.h5' extension")

    analyzer = VideoAnalysis(config)
    calibrator.update_config(config_file)
    scoretype = config_file['scoretype']
    if scoretype == 'single':
        analyzer.analyze(example_behav_file, show_plots = True, continue_analysis = False)
        exit()
    else:
        ##re-load updated config file
        config, rootdir, _, file_list = read_config(config_file) 
        #initialize the program with analysis parameters
        for analyze_file in file_list:
            if analyze_file.endswith('.h5'):
                analyzer.analyze(os.path.join(rootdir, analyze_file), show_plots = False, continue_analysis = True)
        #opportunity to analyze new list with same parameters
        while True:
            keep_going=input("Summary file is complete, would you like to add another folder? (yes/no)  ")
            if keep_going.lower() in ['yes', 'y']:
                rootdir=str(input('What is the new root directory?: '))
                new_list = os.listdir(rootdir)
                for new_file in new_list:
                    analyzer.analyze(new_file)
                continue
            if keep_going.lower() in ['n', 'no', '']:
                print("Exitting the program")
                break

def read_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)
        rootdir = config["rootdir"]
        file = config["example_behavior_file_h5"]
        file_list = os.listdir(rootdir)
        return config, rootdir, file, file_list
    except FileNotFoundError:
        print(f"Error: Configuration file {config_file} not found.")
        sys.exit(1) 
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
        sys.exit(1)
class VideoAnalysis:

    """
    Neatly organizes the program so that it can be easily integrated into other scripts
    """
    def __init__(self, config):

        self.real_sizecm = config["Video_parameters"]["real_sizecm"]
        self.fps = config["Video_parameters"]["fps"]
        self.project_name = config["project_name"]
        self.mcf = config["Arena_parameters"]["movement_correction_factor"]
        self.center_scale = config["Arena_parameters"]["scale_for_center_zone"]
        self.dataframe = None
        self.Frame = None
        self.labels = None
        self.trackBP = config["Arena_parameters"]["bodyparts_to_analyze"]
        self.arena = None
        self.arena1 = None
        self.center_zone = None
        self.border = None
        self.border1 = None
        self.wall = None
        self.length_s = config["Video_parameters"]["length_s"]
        self.startBP = config["Arena_parameters"]["bodyparts_to_determine_trial_start"]
        self.scorer = None
        self.scoretype = config["scoretype"]
        self.likelihood = config["Arena_parameters"]["plikelihood"]
        self.Videos = []; self.Distance=[]; self.Velocity=[]; self.Time_in_Center=[]
        self.Time_in_Border=[]; self.Time_in_Arena=[]; self.BorderEntry=[]; 
        self.Time_immobile = []; self.Time_in_motion = []; self.CenterEntry=[]


    def load_data(self, file):
        """
        loads each file in the folder and extracts the dataframe headers and body part labels"""
        dataframe = pd.read_hdf(file)
        dataframe.head(); scorer = dataframe.columns.get_level_values(0)[0]
        self.dataframe = dataframe[scorer]
        print(self.dataframe['nose']['x'])
        self.labels = self.markers(self.dataframe)

    def get_timestamp(self):
        """
        Grabs current timestamp
        """
        dt = str(datetime.now()).replace(' ', '_').replace(':', '-').split('.')[0]
        return dt
    
    def video_name(self):
        """
        Extracts video name for final dataframe output
        """
        return self.file.replace('DLC', '#').split('#')[0]
    
    def markers(self, x):
        """
        Identifies the bodypart markers that are in each video DataFrame
        """
        bodyparts = [*set(list({i[0] for i in x.columns.values}))]
        return bodyparts
    
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

 
        #shortened:
        # #outter arena
        # arena = np.array([[TLx, TLy], [TRx, TRy], [BRx, BRy], [BLx, BLy]])

        # #center zone
        # difference = self.scale_arena(arena, self.center_scale)
        # grow = self.grow_arena(arena)
        # arena_polygon = Polygon(arena)
        # center_zone = Polygon(arena_polygon.buffer(-difference))
        # arena1 = Polygon(arena_polygon.buffer(grow))
        
        # #establish walls
        # wall = Polygon([Point(p) for p in np.concatenate((center_zone.exterior.coords, arena1.exterior.coords))])
        # border1 = Polygon([Point(p) for p in np.concatenate((center_zone.exterior.coords, arena_polygon.exterior.coords))])
        # border = Polygon([Point(p) for p in np.concatenate((center_zone.exterior.coords, arena_polygon.exterior.coords))])

    def track_labels(self):

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
    
    def start_bp(self):
        print(self.labels)
        k = 0
        self.startBP = []
        
        while k < 3:
            sBP = str(input("What bodyparts (3) should be in the arena at the start of tracking? If only one, enter 3 times "))
            if str(sBP) in self.labels:
                k += 1
                self.startBP.append(str(sBP))
            else:
                print("That bodypart is not in this file")
                continue
        
        print(self.startBP)

    def arena_labels(self):

        print(self.labels)
        k = 0

        while k < 4:
            sBP = str(input("What are the arena labels (4 e.g.)?: "))
            if str(sBP) in self.labels:
                k += 1
                self.arena_labels.append(str(sBP))
            else:
                print("That bodypart is not in this file")
                continue
        
        print(self.arena_labels)

    def crop_video(self):
        """
        crops video based on user input
        
        """

        k=0
        crop=[]
        standard1=self.dataframe[self.startBP[0]]['likelihood'].values
        standard2=self.Dataframe[self.startBP[1]]['likelihood'].values
        standard3=self.Dataframe[self.startBP[2]]['likelihood'].values
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
        crop.append((start+((self.length_s)*self.fps)))
        end=crop[1]
        return start, end;

    def set_up_analysis(self, file):
        self.load_data(file)
        
        # Set trackBP if not already set
        if not self.trackBP:
            self.trackBP = self.track_labels()

        if not self.startBP:
            self.startBP = self.start_bp()

        # Determine score type
        while not self.scoretype:
            scoretype = input("Are you analyzing a single (type 'single') or multiple (type 'batch' or 'multiple') videos?: ").lower()
            if scoretype in ['single', 'batch', 'multiple']:
                self.scoretype = scoretype
            else:
                print("Invalid input. Please enter 'single', 'batch', or 'multiple'.")
                continue
        
        while not self.real_sizecm:
            real_sizecm = input("What is the size of the arena (hypotenuse in centimeters e.g. 65): ")
            if isinstance(real_sizecm, int) and real_sizecm > 30:
                self.real_sizecm = real_sizecm
        
        while not self.fps:
            fps = input("What is the frame rate of the videos to be analyzed?: ")
            if isinstance(fps, int) and real_sizecm > 5:
                self.fps = fps
        # Set center scale
        while not self.center_scale:
            try:
                scale = float(input("Set scale for center zone between 0 and 1 (0.28 is most similar to Ethovision): "))
                if 0 < scale < 1:
                    self.center_scale = scale
                else:
                    print('Scale is outside of range. Please enter a value between 0 and 1.')
            except ValueError:
                print("Invalid input. Please enter a numeric value between 0 and 1.")

        # Create and validate arena
        while True:
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
            
            good = input("Is this arena okay? yes(or ENTER key) or no: ").lower()
            plt.close('all')
            if good in ['yes', 'y', '']:
                break
            elif good in ['no', 'n']:
                continue

    # Crop videos if needed
        while True:
            self.crop = input("Do you want to crop your videos? (yes or no); Start time will be determined by a certain likelihood of the mouse nose: ").lower()
            if self.crop in ['yes', 'y']:
                print(self.startBP)
                while True:
                    try:
                        length_s = int(input("How much time do you wish to analyze in seconds?: "))
                        if length_s > 0:
                            self.length_s = length_s
                            break
                        else:
                            print("Value needs to be greater than 0.")
                    except ValueError:
                        print("Invalid input. Please enter a positive integer.")
                break
            elif self.crop == 'no':
                self.length_s == None
                break
            else:
                print("Invalid input. Please enter 'yes' or 'no'.")

        # Determine likelihood of accuracy
        if not self.likelihood:
            while True:
                try:
                    likelihood = float(input('What is the likelihood cutoff e.g. >0.8. (Input number <1): '))
                    if 0 < likelihood < 1:
                        self.likelihood = likelihood
                        break
                    else:
                        print("Likelihood must be a number between 0 and 1.")
                except ValueError:
                    print("Invalid input. Please enter a number between 0 and 1.")

        # Crop video and analyze
        self.analyze(file, show_plots=True, continue_analysis=False)



    def cleanup(self, x, y, arena1):
        coords = list(zip(x, y))
        cleaned_coords = [(a, b) if Point(a, b).within(arena1) else (np.nan, np.nan) for a, b in coords]
        x_clean, y_clean = zip(*cleaned_coords)
        x_series = pd.Series(x_clean).interpolate(method='linear')
        y_series = pd.Series(y_clean).interpolate(method='linear')
        return x_series, y_series

    def fill_in(self, x, y, likelihood, bp):
        plikelihood = self.dataframe[bp]['likelihood'].values
        coords = list(zip(x, y))
        filled_coords = [(a, b) if lik >= likelihood else (np.nan, np.nan) for (a, b), lik in zip(coords, plikelihood)]
        x_filled, y_filled = zip(*filled_coords)
        x_series = pd.Series(x_filled).interpolate(method='linear')
        y_series = pd.Series(y_filled).interpolate(method='linear')
        return x_series, y_series

    def calculate_distance(self, x, y, arena):
        arenax, arenay = arena.exterior.xy
        dist_arena = float(math.sqrt((arenax[2] - arenax[0]) ** 2 + (arenay[2] - arenay[0]) ** 2))
        real_factor = self.real_sizecm / dist_arena
        coords = list(zip(x[:-1], y[:-1], x[1:], y[1:]))
        distances = [float(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)) * real_factor * self.mcf for x1, y1, x2, y2 in coords]
        return sum(distances)

    def calculate_velocity(self, distance, length_s):
        return distance / length_s

    def calculate_movement(self, x, y):
        coords = list(zip(x[:-1], y[:-1], x[1:], y[1:]))
        move_count = sum(1 for x1, y1, x2, y2 in coords if float(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)) > 0.8)
        no_move_count = len(coords) - move_count
        return move_count / self.fps, no_move_count / self.fps

    def time_in_zones(self, x, y):
        coords = list(zip(x, y))
        center_time = sum(1 for a, b in coords if Point(a, b).within(self.center_zone)) / self.fps
        arena_time = sum(1 for a, b in coords if Point(a, b).within(self.arena1)) / self.fps
        border_time = sum(1 for a, b in coords if Point(a, b).within(self.border)) / self.fps
        return arena_time, center_time, border_time

    def frequency_crossing(self, x, y):
        coords = list(zip(x[:-1], y[:-1], x[1:], y[1:]))
        border_entry = sum(1 for x1, y1, x2, y2 in coords if Point(x1, y1).within(self.center_zone) and Point(x2, y2).within(self.border))
        center_entry = sum(1 for x1, y1, x2, y2 in coords if Point(x1, y1).within(self.border) and Point(x2, y2).within(self.center_zone))
        return border_entry, center_entry

  # Function to load and update the config
    def update_config(self, config_file, updated_values=None):


        yaml = YAML()
        yaml.preserve_quotes = True  # Optional: Preserves quotes around values

        # Load the YAML file with comments
        with open(config_file, 'r') as file:
            config_data = yaml.load(file)
        
        config_data['scoretype'] = self.scoretpye

        #Video_parameters:
        config_data['Video_parameters']['fps'] = self.fps,
        config_data['Video_parameters']['real_sizecm'] = self.real_sizecm
        config_data['Video_parameters']['length_s'] = self.real_sizecm

        #Arena Parameters:
        config_data['Arena_parameters']['bodyparts_to_analyze'] = self.startBP
        config_data['Arena_parameters']['bodyparts_to_determine_trial_start'] = self.trackBP
        config_data['Arena_parameters']['movement_correction_factor'] = self.mcf
        config_data['Arena_parameters']['plikelihood'] = self.likelihood
        config_data['Arena_parameters']['scale_for_center_zone'] = self.center_scale

        # Write the updated data back, preserving comments
        with open(config_file, 'w') as file:
            yaml.dump(config_data, file)

    def analyze(self, file, show_plots = True, continue_analysis = True):
        
        self.file = file
        self.Videos.append(self.video_name())
        self.load_data(file)
        self.make_arena()

        print(self.dataframe)
        if self.length_s is not None:
            start, end = self.crop_video()
            self.dataframe = self.dataframe[start:end]
        

        for bp in self.trackBP:
            x, y = self.dataframe[bp]['x'].values, self.dataframe[bp]['y'].values
            x, y = self.cleanup(x, y, self.arena1)
            x, y = self.fill_in(x, y, 0.9, bp)
            dist = self.calculate_distance(x, y, self.arena)
            vel = self.calculate_velocity(dist, 300)
            self.Distance.append(dist)
            self.Velocity.append(vel)
                    
            move_time, still_time = self.calculate_movement(x, y)
            self.Time_immobile.append(still_time); self.Time_in_motion.append(move_time)
                    
            arena_time, center_time, border_time = self.time_in_zones(x, y)
            border_entry, center_entry = self.frequency_crossing(x, y)

            self.Time_in_Center.append(center_time)
            self.Time_in_Border.append(border_time)
            self.Time_in_Arena.append(arena_time)

            #Frequency in zones#
            border_entry, center_entry=self.frequency_crossing(x, y)
            self.BorderEntry.append(border_entry)
            self.CenterEntry.append(center_entry)
            plt.plot(x,y, label=f'{bp}')

        if show_plots:
            ##Track BP
            plt.legend()
            print(f"Distance: {tuple(zip(self.trackBP, self.Distance))},\n Velocity: {tuple(zip(self.trackBP,self.Velocity))}, \n Time in Center: {zip(self.trackBP,self.Time_in_Center)} \n, Time in Border: {self.Time_in_Border},\n Time in Arena: {self.Time_in_Arena}, \n Border Entries: {self.BorderEntry},\n Center Entries: {self.CenterEntry}")
            plt.show(block=False)
        else:
            plt.close()


        summary_dict = {'Video': file,
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
        if not continue_analysis:
            if self.scoretype.lower() == 'single':
                video_summary=pd.DataFrame(summary_dict)
                summary=os.path.join(self.rootdir, f'{self.video_name(file)}-{self.get_timestamp()}.xlsx')
                video_summary.to_excel(summary, index=False, header=True)
                sys.exit()
        else:
            headersss=[]
            headersss.append(self.Videos)
            headersss.append(self.trackBP)
            ind = pd.MultiIndex.from_product(headersss, names=["video", "label"])
            ind=pd.MultiIndex.from_product
            video_summary=pd.DataFrame(summary_dict, index = ind)

            summary=os.path.join(self.rootdir, f'{self.project_name}-{self.get_timestamp()}.xlsx')
            video_summary.to_excel(summary, index=True, header=True)

    def clear_variables(self):
        self.length_s = self.fps = self.real_sizecm = self.likelihood = self.startBP = self.center_scale = self.mcf = self.trackBP = None

if __name__ == "__main__":
    main(config_file)
