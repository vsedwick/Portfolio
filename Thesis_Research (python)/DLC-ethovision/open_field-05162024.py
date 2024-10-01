import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon
import os
import matplotlib.pyplot as plt
import sys
import math
from datetime import datetime


def main():
    rootdir = "D:\\Exc Dreadds\\dreadds\\OF\\SAL\\Virgin Female"
    file = "D:\\Exc Dreadds\\dreadds\\OF\\CNO\\Males\\VM14_2023-03-27-222442-0000DLC_resnet50_NOR_OpenfieldJan25shuffle1_50000_filtered.h5"
    fps = 30
    real_sizecm = 65
    project_name = 'Saline1_female'
    mcf = 0.58
    scale = 0.28

    analyzer = VideoAnalysis(rootdir, file, fps, real_sizecm, project_name, mcf, scale)
    analyzer.analyze()

class VideoAnalysis:
    def __init__(self, rootdir, file, fps, real_sizecm, project_name, mcf, scale):
        self.rootdir = rootdir
        self.file = file
        self.fps = fps
        self.real_sizecm = real_sizecm
        self.project_name = project_name
        self.mcf = mcf
        self.scale = scale
        self.dataframe = None
        self.labels = None
        self.trackBP = None
        self.arena = None
        self.arena1 = None
        self.center_zone = None
        self.border = None
        self.border1 = None
        self.wall = None

    def load_data(self):
        self.dataframe = pd.read_hdf(self.file)
        self.dataframe.head()
        self.labels = self.markers(self.dataframe)
    
    def get_timestamp(self):
        dt = str(datetime.now()).replace(' ', '_').replace(':', '-').split('.')[0]
        return dt
    
    def video_name(self):
        return self.file.replace('DLC', '#').split('#')[0]
    
    def markers(self, x):
        bodyparts = list({i[1] for i in x.columns.values})
        return bodyparts
    
    def make_arena(self):
        TLx = self.dataframe['TL']['x'].mean()
        TLy = self.dataframe['TL']['y'].mean()
        TRx = self.dataframe['TR']['x'].mean()
        TRy = self.dataframe['TR']['y'].mean()
        BLx = self.dataframe['BL']['x'].mean()
        BLy = self.dataframe['BL']['y'].mean()
        BRx = self.dataframe['BR']['x'].mean()
        BRy = self.dataframe['BR']['y'].mean()

        arena = np.array([[TLx, TLy], [TRx, TRy], [BRx, BRy], [BLx, BLy]])
        difference = self.scale_arena(arena, self.scale)
        grow = self.grow_arena(arena)
        arena_polygon = Polygon(arena)
        center_zone = Polygon(arena_polygon.buffer(-difference))
        arena1 = Polygon(arena_polygon.buffer(grow))
        
        wall = Polygon([Point(p) for p in np.concatenate((center_zone.exterior.coords, arena1.exterior.coords))])
        border1 = Polygon([Point(p) for p in np.concatenate((center_zone.exterior.coords, arena_polygon.exterior.coords))])
        border = Polygon([Point(p) for p in np.concatenate((center_zone.exterior.coords, arena_polygon.exterior.coords))])

        self.arena, self.arena1, self.center_zone, self.border, self.border1, self.wall = arena_polygon, arena1, center_zone, border, border1, wall

    def track_labels(self):
        trackBP = ['TL', 'nose', 'Neck', 'body_center', 'tail_base']
        return trackBP

    def scale_arena(self, arena, scale):
        xs = [i[0] for i in arena]
        ys = [i[1] for i in arena]
        x_center = 0.5 * (min(xs) + max(xs))
        y_center = 0.5 * (min(ys) + max(ys))
        center = Point(x_center, y_center)
        min_corner = Point(min(xs), min(ys))
        return center.distance(min_corner) * scale

    def grow_arena(self, arena):
        xs = [i[0] for i in arena]
        ys = [i[1] for i in arena]
        x_center = 0.5 * (min(xs) + max(xs))
        y_center = 0.5 * (min(ys) + max(ys))
        center = Point(x_center, y_center)
        min_corner = Point(min(xs), min(ys))
        return center.distance(min_corner) * 0.15

    def crop_video(self, startBP, length_s):
        k = 0
        crop = []
        standards = [self.dataframe[bp]['likelihood'].values for bp in startBP]
        for i, vals in enumerate(zip(*standards)):
            if all(v >= 0.9 for v in vals):
                k += 1
            else:
                k = 0
            if k == 3:
                crop.append(i - 2)
                break
        start = crop[0]
        crop.append(start + length_s * self.fps)
        end = crop[1]
        return start, end

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

    def analyze(self):
        self.load_data()
        self.make_arena()
        self.trackBP = self.track_labels()
        start, end = self.crop_video(self.trackBP, 300)
        self.dataframe = self.dataframe[start:end]
        results = []

        for bp in self.trackBP:
            x, y = self.dataframe[bp]['x'].values, self.dataframe[bp]['y'].values
            x, y = self.cleanup(x, y, self.arena1)
            x, y = self.fill_in(x, y, 0.9, bp)
            distance = self.calculate_distance(x, y, self.arena)
            velocity = self.calculate_velocity(distance, 300)
            move_time, still_time = self.calculate_movement(x, y)
            arena_time, center_time, border_time = self.time_in_zones(x, y)
            border_entry, center_entry = self.frequency_crossing(x, y)

            results.append({
                'Bodypart': bp,
                'Total Distance Traveled (cm)': distance,
                'Velocity (cm/s)': velocity,
                'Center Entries': center_entry,
                'Time in Center': center_time,
                'Border Entries': border_entry,
                'Time in Border': border_time,
                'Time Immobile': still_time,
                'Time Moving': move_time,
                'Time in Arena': arena_time
            })

            plt.plot(x, y, label=f'{bp}')
        plt.legend()
        plt.show()

        summary = pd.DataFrame(results)
        timestamp = self.get_timestamp()
        summary.to_excel(os.path.join(self.rootdir, f'{self.project_name}-{timestamp}.xlsx'), index=False)

if __name__ == "__main__":
    main()
