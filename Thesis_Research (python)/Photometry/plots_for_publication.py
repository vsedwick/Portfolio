import matplotlib.pyplot as plt
import pandas as pd

class PublishPlots:
    def __init__(self, config):
        self.config = config
    def pretty_plot(self, color_list, event_dictionary):

        #ETHOGRAM ALIGNED WITH NORMALIZED TRACE
        plt.figure(figsize = (10,4))

        # plt.align_ylabels()


        plt.plot(self.photo_time,self.trace, color = 'k')
        # etho_fig.label_outer()
        plt.ylim([self.config.trace.min()-(0.30*self.trace.min()), self.trace.max()+(0.30*self.trace.max())])
        plt.ylabel(r'$\Delta$F', fontsize = 20)
        plt.xlabel('Time (s)', fontsize = 20)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)

        x_positions = np.arange(0,self.photo_time.max(),(100*self.photo_fps))

        plt.tight_layout()
        count = 0
        stagger = 0
        add_to_legend = []
        for key in event_dictionary:
            add_to_legend.append(key)
            for k in range(len(event_dictionary[key])):
                if event_dictionary[key][k]==1:
                    a=int(k)
                    b=int(k + self.config.post_s * self.photo_fps)
                    c=int(k-self.config.pre_s * self.photo_fps)
                    plt.axvspan(self.photo_time[a], self.photo_time[b], color=color_list[count],  alpha=0.5, ymax = 0.4-stagger, label = key)
                    plt.legend(fontsize = 20)
            count+=1
            stagger+=0.03
                

        plt.show()

    def plot_ethogram_trace(self, event_dictinonary, Representative_image, color_dictionary):
        # Set up the figure and axis for the plot with gridspec_kw to adjust relative heights
        fig, (ax_trace, ax_ethogram) = plt.subplots(2, 1, figsize=(8, 4), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        plt.subplots_adjust(hspace = 0.1)
        ax_trace.tick_params(axis='y', labelsize=12)

        # Plot the photometry trace with a skinnier line
        ax_trace.plot(self.photo_time, self.trace, 'k', linewidth=0.5)  # 'k' for black line
        ax_trace.axvline(self.photo_time[self.photo_start], color = 'k', linestyle = '--')
        # ax_trace.set_title(f'{subject} - {assay}')

        if self.trace_ == 'Scaled_trace' or self.trace_name =='Normalized trace':
            ax_trace.set_ylabel(r'$\Delta$F', size = 15)
        elif self.trace_name == 'Zscore_trace':
            ax_trace.set_ylabel('zF', size = 15)

        # Set the y-limits to be tighter around the data
        trace_min, trace_max = min(self.trace), max(self.trace)
        ax_trace.set_ylim([trace_min - 0.1 * abs(trace_min), trace_max + 0.1 * abs(trace_max)])

        # Plot the ethogram below the photometry trace
        ethogram_base = -0.1  # Start plotting ethogram events at this y-value
        linelengths = 0.1  # Make the event lines shorter
        count = 0
        legend_patches = []
        for i, behavior in enumerate(sorted(color_dictionary.keys())):
            behavior_data = event_dictinonary.get(behavior, [])
            if any(behavior_data):
                # Extract the times at which the behavior occurs
                behavior_times = self.photo__time[behavior_data]
                ax_ethogram.eventplot(behavior_times, lineoffsets=ethogram_base - count, colors=[color_dictionary[behavior]], linelengths=linelengths)
                # Add a label for the behavior
                # ax_ethogram.text(fp_time[-1], ethogram_base - count, behavior, verticalalignment='center', color=colored_behaviors[behavior], fontsize=8)
                patch = mpatches.Patch(color=color_dictionary[behavior], label=behavior)
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

    def selective_event_dict(self, list_of_events):
        fp_behav_scores = {}
        for i in list_of_events:
            placements = [int(p/self.behav_fps * self.photo_fps) for p,_ in zip(range(len(self.config.behav_raw[i])), self.config.behav_raw[i]) if _ == 1]
            fp_array = []
            for frame in range(len(self.photo_frames)):
                if frame in placements:
                    fp_array.append(1)
                else:
                    fp_array.append(0)
            fp_behav_scores.update({f'{i}': np.array(fp_array, dtype = bool)})
        
        return fp_behav_scores
