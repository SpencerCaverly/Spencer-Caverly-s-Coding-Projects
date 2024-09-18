#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:35:49 2024

@author: spencercaverly
"""

#Sandia National Laboratories High Resolution Photocurrent Mapping Software

import matplotlib.pyplot as plt
import numpy as np
import sys

pc_file = '//snlca/home/smcaver/Desktop/Javier_Data_3/07-24-24_JF-II-231/Difficult_Data_Set/Map_250um_Pt-Cr_Collector_RT_50x_Objective_pc.txt'

r_file = '//snlca/home/smcaver/Desktop/Javier_Data_3/07-24-24_JF-II-231/Difficult_Data_Set/Map_250um_Pt-Cr_Collector_RT_50x_Objective_Copy_raman.txt'

noise_ommision_threshold = 0.275 * (10**(-8))
x_coordinate = 1
y_coordinate = 1
target_integer_wave_number = 520
peak_prominence_detection_threshold = 170
mark_peaks = True
restrict_x = False
xlimits = [0,1750]
ylimits = [0,1000]
calibration_mode = True
end_time = 1889
start_time = 165

# Extract Raman Data From txt File

r_lines = []
total_r_lines = 0
with open(r_file) as file:
    for file_line in file:
        r_lines.append(file_line)
        total_r_lines += 1
file.close()
x_positions = []
y_positions = []
wave_numbers = []
r_intensities = []
i = 0
for line in r_lines:
    if i == 0:
        i += 1
    else:
        x_positions.append(float(line.split()[0]))
        y_positions.append(float(line.split()[1]))
        wave_numbers.append(float(line.split()[2]))
        r_intensities.append(float(line.split()[3]))
        
# Determine the Unique X and Y Positions

unique_y_positions = []
unique_y_positions.append(y_positions[0])
latest_y = y_positions[0]
for y_position in y_positions:
    if y_position != latest_y:
        unique_y_positions.append(y_position)
        latest_y = y_position
unique_x_positions = []
unique_x_positions.append(x_positions[0])
latest_x = x_positions[0]
for x_position in x_positions:
    if x_position == unique_x_positions[0] and np.size(unique_x_positions) > 1:
        break
    if x_position != latest_x:
        unique_x_positions.append(x_position)
        latest_x = x_position

# Extract Photocurrent Data From txt File

pc_lines = []
total_pc_lines = 0
with open(pc_file) as file:
    for file_line in file:
        pc_lines.append(file_line)
        total_pc_lines += 1
file.close()
times = []
currents = []
previous = 1000
tail = True
for line in pc_lines:
    if len(line) >= 15:
        if tail == False:
            times.append(float(line.split()[1]))
            currents.append(float(line.split()[3]))
            if times[-1] >= end_time:
                break
            if line.split()[1].isnumeric() and tail == True and float(line.split()[1]) >= start_time:
                if float(line.split()[3]) > (previous*2):
                    times.append(float(line.split()[1]))
                    currents.append(float(line.split()[3]))
                    tail = False
                previous = float(line.split()[3])
            
# Plot Noisy Photocurrent Data to Determine Proper Noise Ommision Threshold

plt.plot(times, currents)
plt.title("Photocurrent Without Noise Ommision")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.show()

# Filter Out Noise

points_size = np.size(times)
currents_ave = np.mean(currents)
filtered_currents = []
i = 0
for current in currents:
    if current < noise_ommision_threshold and current >= 0:
        filtered_currents.append(current)
    else:
        filtered_currents.append(currents_ave)

# Plot Filtered Photocurrent Data

plt.plot(times, filtered_currents)
plt.title("Filtered Photocurrent")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.ylim(-max(filtered_currents), 2*max(filtered_currents))
plt.show()

# Extract Target Wave Number Indexes and Corresponding Intensity Values

wave_number_indexes = []
k = 0
wave_number_indexes = []
for wave_number in wave_numbers:
    if np.trunc(wave_number) == float(target_integer_wave_number):
        wave_number_indexes.append(k)
    k += 1
if np.size(wave_number_indexes) != (np.size(unique_x_positions) * np.size(unique_y_positions)):
    print("Target integer wave number of " + str(target_integer_wave_number) + " not found in txt file. Either change the number slightly or check txt file for valid inputs.")
    print("Generally, valid input values are within 3 cm^-1 of the desired wave number value.")
    print("Suggested input values: " + str(target_integer_wave_number - 2) + ", " +
str(target_integer_wave_number - 1) + ", " + str(target_integer_wave_number + 1) + ", " +
str(target_integer_wave_number + 2))
    sys.exit()
target_intensities = []
for index in wave_number_indexes:
    target_intensities.append(r_intensities[index])
    
# Extract Photocurrent Values

pc_scans_size = points_size/np.size(unique_y_positions)
point_readings_size = (points_size/np.size(unique_y_positions))/np.size(unique_x_positions)
point_scan = []
local_mins = []
local_maxes = []

# Model 3 - Determining Point Scan Partiton Time Values

previous_current = filtered_currents[0]
partition_time_values = []
temp_array = []
append_temp = False
for (time, filtered_current) in zip(times, filtered_currents):
    if filtered_current <= (0.75 * previous_current):
        append_temp = True
    if append_temp == True:
        temp_array.append(time)
    if filtered_current >= (2 * previous_current) and append_temp == True:
        partition_time_values.append(temp_array[int(np.trunc(np.size(temp_array)/2))] * 10)
        temp_array = []
        append_temp = False
    previous_current = filtered_current
point_spacings = []
for i in range(1, np.size(partition_time_values)):
    point_spacings.append(partition_time_values[i] - partition_time_values[i-1])
average_point_spacing = np.mean(point_spacings)
for i in range(np.size(partition_time_values)-1):
    if partition_time_values[i+1] - partition_time_values[i] <= (0.6 * average_point_spacing):
        partition_time_values[i] = -1
k = 0
for i in range(np.size(partition_time_values)):
    i -= k
    if partition_time_values[i] <= 0:
        del(partition_time_values[i])
        k += 1
if np.size(partition_time_values) < (np.size(unique_x_positions) * np.size(unique_y_positions)):
    point_spacings = []
    for i in range(1, np.size(partition_time_values)):
        point_spacings.append(partition_time_values[i] - partition_time_values[i-1])
    for i in range(1, np.size(partition_time_values)):
        if (partition_time_values[i] - partition_time_values[i-1]) == max(point_spacings):
            partition_time_values.insert(partition_time_values.index(partition_time_values[i]), (((partition_time_values[i] - partition_time_values[i-1])/2) + partition_time_values[i-1]))
            del(point_spacings[point_spacings.index(max(point_spacings))])
            if np.size(partition_time_values) == (np.size(unique_x_positions) * np.size(unique_y_positions)):
                break
i = start_time * 10 # pc points
j = 0 # partition index
for filtered_current in filtered_currents:
    if i <= partition_time_values[j]:
        point_scan.append(filtered_current)
    if i > partition_time_values[j]:
        local_maxes.append(max(point_scan))
        local_mins.append(min(point_scan))
        point_scan = [filtered_current]
        j += 1
    i += 1
    if j == np.size(partition_time_values) or i >= end_time * 10:
        break
dark_current = max(local_mins)
photocurrents = []
row = []
j = 0
i = 0

for local_max in local_maxes:
    if i == 0:
        photocurrents.append(local_max-dark_current)
        j += 1
    if i >= 1:
        row.append(local_max-dark_current)
        j += 1
    if j == np.size(unique_x_positions):
        if np.size(row) == 0:
            i += 1
            j = 0
    if np.size(row) != 0:
        photocurrents = np.vstack([photocurrents, row])
        row = []
        j = 0
        
# Matricize Raman Data

r_matrix = []
row = []
for i in range(np.size(target_intensities)):
    if i < np.size(unique_x_positions):
        r_matrix.append(target_intensities[i])
    if i/np.size(unique_x_positions) - np.fix(i/np.size(unique_x_positions)) == 0 and i != np.size(unique_x_positions) and i != 0:
        r_matrix = np.vstack([r_matrix,row])
        row = []
    if i == (np.size(target_intensities) - 1):
        row.append(target_intensities[-1])
        r_matrix = np.vstack([r_matrix,row])
    if i >= np.size(unique_x_positions):
        row.append(target_intensities[i])
        
# Determine Raman Intensities and Local Peak Locations

if x_coordinate > np.size(unique_x_positions):
    print("Requested coordinate does not fit within the bounds of the scan.")
    print("Pick a x_coordinate value that is <= " + str(np.size(unique_x_positions)))
    sys.exit()
if y_coordinate > np.size(unique_y_positions):
    print("Requested coordinate does not fit within the bounds of the scan.")
    print("Pick a y_coordinate value that is <= " + str(np.size(unique_y_positions)))
    sys.exit()
raman_intensities = []
i = 0
raman_range = 0
for intensity in r_intensities:
    if y_positions[i] == unique_y_positions[y_coordinate-1] and x_positions[i] == unique_x_positions[x_coordinate-1]:
        raman_intensities.append(intensity)
        raman_range += 1
    i += 1
local_peaks = []
noise_intensity = np.mean(raman_intensities[0:100])
prior_prominent_intensity = 0
local_peak_found = False
for raman_intensity in raman_intensities:
    if raman_intensity < prior_prominent_intensity and local_peak_found == False and prior_prominent_intensity >= (peak_prominence_detection_threshold + noise_intensity):
        local_peaks.append(prior_prominent_intensity)
        local_peak_found = True
    if raman_intensity >= prior_prominent_intensity and local_peak_found == True:
        local_peak_found = False
    prior_prominent_intensity = raman_intensity
local_peak_locations = []
for local_peak in local_peaks:
    index = raman_intensities.index(local_peak)
    local_peak_locations.append(wave_numbers[0:raman_range][index])
    
# Plot Raman Intensity Map

plt.pcolormesh(r_matrix)
plt.colorbar()
plt.title("Raman Intensity Map at " + str(target_integer_wave_number) + " Inverse Centimeters")
plt.plot(x_coordinate - 0.5,y_coordinate - 0.5,"r.")
plt.show()

# Plot Raman Spectra at User Specified Point

if mark_peaks == True:
    for local_peak_location in local_peak_locations:
        if local_peak_location == local_peak_locations[0]:
            plt.plot(local_peak_location*np.ones(2),[0,(ylimits[1]-ylimits[0])*(0.6)],"green", label = "Peak Location")
        else:
            plt.plot(local_peak_location*np.ones(2),[0,(ylimits[1]-ylimits[0])*(0.6)],"green")
plt.plot(wave_numbers[0:raman_range],raman_intensities,"blue")
plt.title("Raman Spectra at point (" + str(x_coordinate) + "," + str(y_coordinate) + ")")
plt.xlabel("Raman Shift (cm^-1)")
plt.ylabel("Intensity")
plt.legend(loc = "upper left")
if restrict_x == True:
    plt.xlim(xlimits[0],xlimits[1])
plt.ylim(ylimits[0],ylimits[1])
plt.show()
print("Peak Locations Detected at " + str(local_peak_locations))

# Plot Photocurrent Intensity Map

if calibration_mode == False:
    plt.pcolormesh(photocurrents)
    plt.colorbar()
    plt.title("Photocurrent Intensity Map (A)")
    plt.show()
    
    
    
# Notes:
# Individual signals with too low of an intensity will have noise that will be near that of the max signal stength or higher, not to exceed user specified noise ommission threshold
# This causes the first column to be noisy and display highly amplified intensity compared to the actual photocurrent value



if calibration_mode == True:
    for i in range(np.size(unique_y_positions)):
        plt.plot(times, filtered_currents, "blue")
        for partition_time_value in partition_time_values:
            plt.plot([partition_time_value/10, partition_time_value/10], [0, max(filtered_currents)], "red")
        plt.xlim((pc_scans_size/10)*i, (pc_scans_size/10)*(i+1))
        plt.title("Graph " + str(i+1))
        plt.show()






