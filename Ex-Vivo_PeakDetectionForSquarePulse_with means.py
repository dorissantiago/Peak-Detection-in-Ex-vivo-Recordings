# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 17:18:28 2025

@author: Doris Santiago
"""

# -*- coding: utf-8 -*-
"""Created on Tue Jun 16 11:32:37 2025@author: Doris Santiago: corrected for gaussian filter"""
"""Created on Wed Sep 25 16:15:58 2024@author: Doris Santiago """
"""Modified on Mon Feb  4 16:15:58 2025 @author: Doris Santiago"""
"""Modified on Wed Mar 26 15:42:58  2025 @author: Doris Santiago"""

# FUNCTIONAL CODE FOR PEAK DETECTION


import pyabf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os

from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d 


# =============================================================================
# #test files
# fn = r'C:\Users\Doris Santiago\Documents\Experiments\Ex-vivo\EphysAnalysis\DataForAnalysis\RawData\gabazine\2024-11-13\B7_096_MSN7_SqPulse_GFP05_0000.abf'
# # fn = r'C:\Users\Doris Santiago\Documents\Experiments\Ex-vivo\EphysAnalysis\DataForAnalysis\RawData\gabazine\2024-11-29\B8_265_MSN44_IVCurveEnd_0002.abf'
# # fn = r'C:\Users\Doris Santiago\Documents\Experiments\Ex-vivo\EphysAnalysis\DataForAnalysis\RawData\gabazine\2024-11-13\B7_096_MSN7_SqPulse_GFP05_0000.abf'
# # fn= r'C:/Users/Doris Santiago/Documents/Experiments/Ex-vivo/EphysAnalysis/DataForAnalysis/RawData/gabazine/2025-02-20/B11_051_siteS1_MSN9_GFP35_STP_0001.abf'
# fn = r'C:/Users/Doris Santiago/Documents/Experiments/Ex-vivo/EphysAnalysis/DataForAnalysis/RawData/Gabazine/2025-03-04/B12_057_siteM1_MSN18_GFP60_0000.abf'

# =============================================================================
# =============================================================================
# add the file to analyse here



fn= r'add file path here '


# # # # =============================================================================
abf = pyabf.ABF(fn)
fs = abf.sampleRate
labelname =abf.abfID
# print(labelname)

#Inititate the plot 
fig, ax = plt.subplots(1, figsize=(7, 4), dpi=200)
# =============================================================================
# DEFINE BASELINE AND EXTRACT EPOCHS
# =============================================================================
SD = 5
prominence_value = 0.5
 #to be used in two loops
sigma = 0.1


df_peaks = pd.DataFrame()
for sweepNumber in abf.sweepList:
    abf.setSweep(sweepNumber)
    sigma = sigma # Change this if needed
    filtered_sweepY = gaussian_filter1d(abf.sweepY, sigma)
    time = abf.sweepX  # time for the current sweep
    baseline_duration_ms = 100
    sampling_rate = abf.sampleRate
    digital_output = abf.sweepD(0)
    
    try:
        baseline_samples = int((baseline_duration_ms / 1000) * sampling_rate)
        epoch_start = abf.sweepEpochs.p1s[2]
        # epoch_start = abf.sweepEpochs.p1s[4]
        baseline_start_index = int(epoch_start) - baseline_samples
        baseline_end_index = epoch_start
        
        # Calculate RMP as mean of voltage in baseline
        baseline_voltage = np.mean(filtered_sweepY[baseline_start_index:baseline_end_index])
        bw = range(baseline_start_index, baseline_end_index + 1)
        baseline_voltage = np.mean(filtered_sweepY[bw]) if bw else np.nan  # Mean current over baseline period

        epoch_start_points = abf.sweepEpochs.p1s
        epoch_end_points = abf.sweepEpochs.p2s
        epoch_start_points = [p for p in epoch_start_points if p < len(digital_output)]
        epoch_end_points = [p for p in epoch_end_points if p < len(digital_output)]
        
        # Iterate through each epoch to create the DataFrame
        data = []
        for i, (start, end) in enumerate(zip(epoch_start_points, epoch_end_points)):           
            epoch_output = digital_output[start:end]
            if all(val in [0, 1] for val in epoch_output):
                data.append({
                    'EpochStart': start,
                    'EpochEnd': end,
                    'DO_Start': digital_output[start],
                    'DO_End': digital_output[end],
                    'epoch_time': end - start,
                    'epochType': abf.sweepEpochs.types[i],
                })

        # Create DataFrame from the data list
        df_epochs = pd.DataFrame(data)
        # print(df_epochs)
        
        extracted_df_epochs = df_epochs[df_epochs['DO_Start'] == 1]      # 1 means light ON . 0 means light OFF
        LightON = extracted_df_epochs[['EpochStart', 'EpochEnd', 'DO_Start','DO_End','epoch_time']]
        print(LightON)

        stimulation_duration_seconds = extracted_df_epochs['epoch_time']/sampling_rate
               
    except Exception as e:
        print(f"An error occurred: {e}")
    
    # recovery_time = df_epochs['epoch_time'].iloc[5]
    
# =============================================================================
# FIND PEAKS
# =============================================================================


peak_data_list =[]
for sweep_num in range(abf.sweepCount):
    abf.setSweep(sweep_num)
    sigma = sigma # Change this if needed
    filtered_sweepY = gaussian_filter1d(abf.sweepY, sigma)
    time = abf.sweepX   # time for the current sweep
    baseline_duration_ms = 100
    sampling_rate = abf.sampleRate
    
    try:
        baseline_samples = int((baseline_duration_ms / 1000) * sampling_rate)
        # # epoch_start = abf.sweepEpochs.p1s[2]
        # epoch_start = abf.sweepEpochs.p1s[4]
        epoch_start = extracted_df_epochs["EpochStart"].iloc[0],
        baseline_start_index = int(epoch_start) - baseline_samples
        baseline_end_index = epoch_start
        
        # Calculate RMP as mean of voltage in baseline
        baseline_voltage = np.mean(filtered_sweepY[baseline_start_index:baseline_end_index])
        bw = range(baseline_start_index, baseline_end_index)
        baseline_voltage = np.mean(filtered_sweepY[bw]) if bw else np.nan  # Mean current over baseline period
    except:
       pass
        
    peak_threshold = baseline_voltage  + SD *(np.std(filtered_sweepY[baseline_start_index:baseline_end_index]))   # Adjust threshold based on RMP
    peak_prominence = prominence_value # Adjust prominence (CAREFUL !!! lower the prominence for smaller peaks)   Standard prominence is 5
    min_distance = 20  # Minimum distance between peaks (adjust based on your data)
    
    # Find peaks with `find_peaks` (no need to request widths here)
    peaks, properties = find_peaks(
        filtered_sweepY,
        prominence=peak_prominence,
        distance=min_distance
    )

    # Initialize a list to store peak currents for the current sweep
    peak_currents = []

    # Extract peak properties
    for i, peak_idx in enumerate(peaks):
        peak_time = abf.sweepX[peak_idx]
        peak_voltage = filtered_sweepY[peak_idx]  # Assuming abf.sweepC holds the voltage trace
        
        ###  Find Peak Threshold  ###
        #Calculate the first derivative
        dVdt = np.gradient(filtered_sweepY, abf.sweepX)  # dV/dt calculation  
        baseline_dVdt = dVdt[baseline_start_index:baseline_end_index]  # Baseline period
        threshold_dVdt = np.mean(baseline_dVdt) + SD * np.std(baseline_dVdt)  # Threshold for significant depolarization
        
    
        # Find the first point where dV/dt exceeds threshold
        threshold_idx = np.where(dVdt > threshold_dVdt)[0][0]
        peak_threshold_voltage =filtered_sweepY[threshold_idx]
        peak_threshold_time = abf.sweepX[threshold_idx]
            
        ## Peak Amplitude ##  (Peak - Baseline Voltage)
        peak_amplitude = abs(baseline_voltage - peak_voltage)     
        
        ###  Peak Width  ###
        half_max_voltage = baseline_voltage + (peak_amplitude / 2)                              # Compute Half-Max Voltage
        crossings = np.where(np.diff(np.sign(filtered_sweepY - half_max_voltage)))[0]           # Find where the voltage crosses the half-max voltage
        print(f"Crossings indices: {crossings}")
        
        width_start_idx = crossings[0]
        width_end_idx = crossings[1]
        
        width_s = (abf.sweepX[width_end_idx]) - (abf.sweepX[width_start_idx])
        width_ms = width_s * 1000           #Convert s to ms

        ### Find UpSlope and DownSlope  ###
        # Upslope is calculated as the slope between the rising or depolarization width to the peak voltage
        us_V_start = filtered_sweepY[width_start_idx]
        us_V_peak = filtered_sweepY[peak_idx]               #this is the same as the peak voltage        
        us_T_start = abf.sweepX[width_start_idx]
        us_T_peak = abf.sweepX[peak_idx]
        
        peak_upslope = (us_V_peak - us_V_start) / ((us_T_peak - us_T_start) * 1000)  # in mV/ms # Convert to mV/ms
        # print(f"Upslope: {peak_upslope:.2f} mV/ms")
        
        #Downslope is calculated as the slope between the peak and falling phase or repolarization until half width
        ds_V_peak = filtered_sweepY[peak_idx]
        ds_V_end = filtered_sweepY[width_end_idx]  # Voltage at the point where width ends
        ds_T_peak = abf.sweepX[peak_idx]
        ds_T_end = abf.sweepX[width_end_idx]
        
        peak_downslope = (ds_V_peak - ds_V_end) / ((ds_T_peak - ds_T_end) * 1000)          # in mV/ms # Convert to mV/ms
        
        # Inter-spike interval (time between this peak and the next peak)
        if i + 1 < len(peaks):
            next_peak_idx = peaks[i + 1]
            inter_spike_interval = abf.sweepX[next_peak_idx] - peak_time
        else:
            inter_spike_interval = np.nan  # No next peak, set as NaN
            
        # Append the current peak's current value to the list for this sweep
        clean_path = fn.replace('\\', '/')
        
        peak_currents.append(peak_currents)

        # Append to list
        peak_data_list.append({
            "col1":None,
            "col2": None,
            "FileID": labelname,
            "Sweep_Number": sweep_num + 1,
            "Peak_Number": i + 1,
            "Start_of_Pulse":(extracted_df_epochs["EpochStart"].iloc[0])/abf.sampleRate,
            # "Start_of_Pulse": abf.sweepEpochs.p1s[4]/fs,
            "Peak_Time(s)": peak_time,
            "Peak_Voltage(mV)": peak_voltage,
            "Peak_Amplitude (mV)" : peak_amplitude,
            "Peak_Width(ms)":width_ms,
            "Peak_Onset(s)":peak_threshold_time,
            "Peak_Threshold_Voltage(mV)":peak_threshold_voltage,
            "Peak_Upslope (mV/ms)":peak_upslope,
            "Peak_Downslope(mV/ms)": peak_downslope,
            "BaselineVoltage": baseline_voltage,
            "FilePath": clean_path,
            "Comment": None
            # "Interspike_interval": inter_spike_interval
            })
               
    df_peaks = pd.DataFrame(peak_data_list)
    # print(df_peaks)
            
    #PLOT THE DATA
    ax.set_title(f"Sweep {sweep_num}, Baseline Voltage: {baseline_voltage:.2f} mV, Peak Voltage: {peak_voltage:.2f} mV, SD = {SD}", fontsize = 10)
    
    ax.plot(abf.sweepX, filtered_sweepY, color='#1f77b4')              #Plot Time vs filtered sweepY
    ax.plot(time[peaks], filtered_sweepY[peaks], "go", label="Peaks")  # Mark detected peaks
    
    # Mark epochs with vertical dashed lines
    for p1 in abf.sweepEpochs.p1s:
        ax.axvline(time[p1], color='k', ls='--', alpha=0.5)
     
    #Draw Baseline
    ax.axhline(y= baseline_voltage, color='black', linestyle='--', label='Mean Voltage (baseline_voltage)')
    ax.plot(time[threshold_idx], peak_threshold_voltage, 'o', color='gold', label='Peak Threshold', markersize=8)
    
    ax.set_xlabel(abf.sweepLabelX)
    ax.set_ylabel(abf.sweepLabelY)
    ax.set_ylim(-90, -55)
    ax.set_xlim(0.3, 0.8)
    # ax.set_xlim(0.25, 0.65)
    plt.xticks( fontsize=10)


#end the loop with ending the plot
fig.suptitle(f"{labelname}, Protocol: {abf.protocol}, GaussianFilter(σ):{sigma}", fontsize=10, )
plt.legend(fontsize = 5, loc = "upper right")
plt.tight_layout()
plt.show()

save_name1 = f"{labelname}_allSweeps_{abf.protocol}.png"
save_path1 = os.path.join(r'D:\EphysAnalysis\2025\03_SquarePulseAnalysis\3B\allImages', save_name1)   #REPLACE TO THE PATH WHERE YOU WANT TO SAVE THE IMAGES
plt.savefig(save_path1)
print(f"Image saved to {save_path1}")

# =============================================================================
# # Compute mean trace across all sweeps
# =============================================================================
all_sweeps = []  # List to store all sweep data

for sweep_num in range(abf.sweepCount):
    abf.setSweep(sweep_num)
    sigma = sigma # Change if needed
    filtered_sweepY = gaussian_filter1d(abf.sweepY, sigma)
    all_sweeps.append(filtered_sweepY)
all_sweeps = np.array(all_sweeps)
mean_trace = np.mean(all_sweeps, axis=0)  # Mean trace of all sweeps

#mean baseline voltage
mean_bw = range(baseline_start_index, baseline_end_index)
mean_baseline_voltage = np.mean(mean_trace[mean_bw]) if mean_bw else np.nan
 
peak_prominence = prominence_value  # Adjust prominence as needed
min_distance = 20  # Minimum distance between peaks

# Find peaks in mean trace
peaks, properties = find_peaks(
    mean_trace,
    prominence=peak_prominence,
    distance=min_distance
)
# Initialize list for mean peak properties
mean_peak_data = []  # List to store data for each peak

for i, mean_peak_idx in enumerate(peaks):
    mean_peak_time = abf.sweepX[mean_peak_idx]
    mean_peak_voltage = mean_trace[mean_peak_idx]
    
    ###  Find Peak Threshold  ###
    # Calculate the first derivative of the mean trace
    mean_dVdt = np.gradient(mean_trace, abf.sweepX)  # dV/dt calculation
    mean_baseline_dVdt = mean_dVdt[baseline_start_index:baseline_end_index]  # Baseline period
    mean_threshold_dVdt = np.mean(mean_baseline_dVdt) + SD * np.std(mean_baseline_dVdt)  # Threshold for significant depolarization

    # Find the first point where dV/dt exceeds threshold
    mean_threshold_idx = np.where(mean_dVdt > mean_threshold_dVdt)[0][0]
    mean_peak_threshold_voltage = mean_trace[mean_threshold_idx]
    mean_peak_threshold_time = abf.sweepX[mean_threshold_idx]
    
    # Peak Amplitude : Peak - Baseline Voltage
    mean_peak_amplitude = abs(mean_baseline_voltage - mean_peak_voltage)  # Absolute difference from baseline
    
    ### Peak Width ###
    mean_half_max_voltage = mean_baseline_voltage + (mean_peak_amplitude / 2)  # Compute Half-Max Voltage
    mean_crossings = np.where(np.diff(np.sign(mean_trace - mean_half_max_voltage)))[0]  # Find where the voltage crosses the half-max voltage
    
    if len(mean_crossings) >= 2:
        mean_width_start_idx = mean_crossings[0]
        mean_width_end_idx = mean_crossings[1]
    
        mean_width_s = (abf.sweepX[mean_width_end_idx]) - (abf.sweepX[mean_width_start_idx])
        mean_width_ms = mean_width_s * 1000  # Convert from seconds to milliseconds
    else:
        mean_width_ms = np.nan  # If no valid crossing is found
    
    ### Find UpSlope and DownSlope ###
    # Upslope: Slope between rising phase and peak voltage
    mean_us_V_start = mean_trace[mean_width_start_idx]
    mean_us_V_peak = mean_trace[mean_peak_idx]  # Same as peak voltage        
    mean_us_T_start = abf.sweepX[mean_width_start_idx]
    mean_us_T_peak = abf.sweepX[mean_peak_idx]
    
    mean_peak_upslope = (mean_us_V_peak - mean_us_V_start) / ((mean_us_T_peak - mean_us_T_start) * 1000)  # in mV/ms
    
    # Downslope: Slope between peak and falling phase to half width
    mean_ds_V_peak = mean_trace[mean_peak_idx]
    mean_ds_V_end = mean_trace[mean_width_end_idx]  # Voltage at the point where width ends
    mean_ds_T_peak = abf.sweepX[mean_peak_idx]
    mean_ds_T_end = abf.sweepX[mean_width_end_idx]
    
    mean_peak_downslope = (mean_ds_V_peak - mean_ds_V_end) / ((mean_ds_T_peak - mean_ds_T_end) * 1000)  # in mV/ms
    
    # Inter-spike interval (time between this peak and the next peak)
    if i + 1 < len(peaks):
        mean_next_peak_idx = peaks[i + 1]
        mean_inter_spike_interval = abf.sweepX[mean_next_peak_idx] - mean_peak_time
    else:
        mean_inter_spike_interval = np.nan  # No next peak, set as NaN
    
    # Append the mean peak data for the current peak
    mean_peak_data.append({
        "col1":None,
        "col2": None,
        "FileID": labelname,
        "Sweep_Number": "mean",
        "Peak_Number": i + 1,
        "Start_of_Pulse": (extracted_df_epochs["EpochStart"].iloc[0])/abf.sampleRate,
        # "Start_of_Pulse": abf.sweepEpochs.p1s[2] / abf.sampleRate,
        # "Start_of_Pulse": abf.sweepEpochs.p1s[4] / abf.sampleRate,
        "Peak_Time(s)": mean_peak_time,
        "Peak_Voltage(mV)": mean_peak_voltage,
        "Peak_Amplitude (mV)": mean_peak_amplitude,
        "Peak_Width(ms)": mean_width_ms,
        "Peak_Onset(s)": mean_peak_threshold_time,
        "Peak_Threshold_Voltage(mV)": mean_peak_threshold_voltage,
        "Peak_Upslope (mV/ms)": mean_peak_upslope,
        "Peak_Downslope(mV/ms)": mean_peak_downslope,
        "BaselineVoltage": mean_baseline_voltage,
        # "Interspike_interval": mean_inter_spike_interval,
        "FilePath": clean_path,
        "Comment": None})

# Convert to DataFrame
df_mean_peaks = pd.DataFrame(mean_peak_data)
# print(df_mean_peaks)

# Plot the mean trace with detected peaks
plt.figure(figsize=(13, 5))
plt.plot(abf.sweepX, mean_trace, label="Mean Trace", color="blue")
plt.plot(abf.sweepX[peaks], mean_trace[peaks], "ro", label="Detected Peaks")  # Mark peaks

# Plot baseline voltage
plt.axhline(mean_baseline_voltage
            , color="black", linestyle="--", label="Baseline Voltage")
# plt.axhline(baseline_voltage, color="black", linestyle="--", label="Baseline Voltage")


# Plot Peak Threshold
plt.plot(abf.sweepX[mean_threshold_idx], mean_peak_threshold_voltage, 'o', color='gold', label='Peak Threshold', markersize=8)

# Mark epoch start points with vertical dashed lines
for p1 in abf.sweepEpochs.p1s:
    plt.axvline(abf.sweepX[p1], color='k', ls='--', alpha=0.5)

plt.title(f"Mean Trace of {labelname} in {abf.protocol},Mean Amplitude size = {mean_peak_amplitude:.2f}mV, Mean Baseline voltage = {mean_baseline_voltage:.2f}mV", fontsize= 11)
plt.xlabel("Time (s)")
plt.ylabel("Voltage (mV)")
# plt.xlim(0.3, 0.6)
# plt.xlim(0.25, 0.65)
# plt.ylim(-90,-50)
plt.legend(fontsize = 8)
plt.tight_layout()
plt.show()

# # =============================================================================
# # #PLOT TO TEST THE GAUSIAN FILTER
# # =============================================================================
# Create subplots

# # Create a single plot
# fig, ax = plt.subplots(figsize=(12, 6))

# ax.set_title(f'Raw and Gaussian Filtered Traces (σ={sigma})')
# ax.plot(abf.sweepX, abf.sweepY, color='b', label='Raw Trace')  # Raw trace
# ax.plot(abf.sweepX, filtered_sweepY, color='r', label=f'Gaussian Filter (σ={sigma})')  # Filtered trace
# for p1 in abf.sweepEpochs.p1s:
#     ax.axvline(abf.sweepX[p1], color='k', ls='--', alpha=0.5)
# ax.legend()

# # Adjust layout for a cleaner look
# plt.tight_layout()

# # Show the plot
# plt.show()

#============================================================================
# SAVING EVERTHING
# =============================================================================
#REPLACE TO THE PATH WHERE YOU WANT TO SAVE THE ALL THE DIFERENT FILES

result_folder = r'D:\EphysAnalysis\2025\03_SquarePulseAnalysis\3B\allExcel'
protocol_name = abf.protocol  # Replace with the actual protocol name from your abf object

save_name = f"{labelname}_{protocol_name}.png"
result_folder_images = r'D:\EphysAnalysis\2025\03_SquarePulseAnalysis\3B\allImages'
save_path = os.path.join(result_folder_images, save_name)
plt.savefig(save_path)
print(f"Image saved to {save_path}")

output_path3 = os.path.join(result_folder, f"{labelname}_df_peaks_SquarePulse2ms.csv")
df_peaks.to_csv(output_path3, index=False)
print(f"Peak data saved to {output_path3}")

output_path4 = os.path.join(result_folder, f"{labelname}_df_epochs_SquarePulse2ms.csv")
df_epochs.to_csv(output_path4, index=False)
print(f"Epoch data saved to {output_path4}")

output_path5 = os.path.join(result_folder, f"{labelname}_df_mean_peaks_SquarePulse2ms.csv")
df_mean_peaks.to_csv(output_path4, index=False)
print(f"Peak data saved to {output_path3}")