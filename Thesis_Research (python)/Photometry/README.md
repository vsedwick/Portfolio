# Fiber Photometry Analysis Scripts

This collection of scripts began my coding journey as my lab did not have a sufficient analysis pipeline for the bulk of data we needed to analyze. These were some of the earliest and most complex scripts I have written. 

In present day, I am actively cleaning the scripts to improve readibility and organization. I am addressing issues such as verbose functions and redundancies and adding 
descriptive docstrings. The goal is to publish this collection as an open source analysis program to aid researchers in their work.

The benefit to this program is that peri-event extraction is automatic significantly reducing user-bias and ensuring uniformity and organization across samples. Additionally, it is capable of batch processing any number of recording trials with minimal, prompted supervision, and pools relevant data from all of the samples into a single csv file.

## Description of Program
### Photocode scripts
The photocode scripts perform the bulk of the data cleaning and data extraction. **"[batch_photocode_vs.py](https://github.com/vsedwick/Portfolio/blob/master/Thesis_Research%20(python)/Photometry/batch_photocode_v2.py)"** is the most recent version and is still in the process of being updated whereas "combined-photocode-05052024.py" and other variations are the versions I worked on and utilized upwards of 5 years in my lab. All extractions are nestled neatly into a folder called "Behaviors" and are organized by the peri-event they belong to. Approximately 8 parameters are extracted per sample, per event. These values can be accessed or re-analyzed at anytime without having to re-run this 'photocode'.

### Summary Scripts
The summary script (summary_values_05042024.py) processes all extracted peri-event parameters from the "Behavior" folder and calculates relevant metrics for each one. These metrics are saved as individual CSV files, categorized by the type of calculation performed (e.g., Δ area under the curve, magnitude change, z-score, etc.). Each CSV includes the sample ID and the calculated value for each observed event. Additionally, averaged z-score plots with shaded standard error of the mean (SEM) are generated for each event, providing an immediate visualization of the overall signal response.

### Statistical Scripts
- VS_Photo_LMM_betweengroups_BATCH.py

### Visualization Scripts
- zscoreFigures_andBins_change colors.py
- pretty_picture_03253034.py

### Utility Scripts
- Compile_wID.py

## Instructions
- To run the batch_photocode, the compatible data can be found on Google Drive:
https://drive.google.com/drive/folders/1JRW4vdPXH9-fDm-3SlwxG2M2b4yX3qcl?usp=sharing

- Set the base folder as the project directory in the configuration file ("config.yaml"). Choose one of the samples to use as the example behavior file.

- After running the photocode, "summary_values_05052024.py" will perform a series of calculations to output plottable values into csv sheets. [This script is currently compatible with the config file, but it will be soon.
