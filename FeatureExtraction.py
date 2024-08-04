# pain profile research
# 2183550
# 25/7/22

## When i wrote this code, only myself and God knew how it worked. And now, only God knows

# input file name as a string (to run a specific file)

# Import the required libraries
from math import floor
import os  # operating system commands
os.chdir("D:\Pain Serum Calcium Analysis files")
import numpy as np  # numerical python functionality
import pandas as pd
import matplotlib.pyplot as plt
import csv
from itertools import zip_longest
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.integrate import trapezoid
from scipy import signal
from scipy.signal import find_peaks_cwt
from scipy.signal import peak_widths
from scipy.integrate import quad

input_main_folder = "Data"
input_growth_factor = "GDNF All"


input_trial_run = ""


input_folder_name = input_main_folder + "/" + input_growth_factor
input_file_name = input_main_folder + "/" + input_growth_factor + "/" + input_trial_run + ".csv"


with open(str(input_folder_name +  '_feature_extraction.csv'), 'w', newline='') as f:              # makes a new file for the extracted data to be placed in
    thewriter = csv.writer(f)
    thewriter.writerow([
        "EC", \
        "Cell Type", \
        "Serum Type", 
        "Sex", \
        "ROI", \
        "MIF", \
        "TMIF", \
        "FDM", \
        "TFDM", \
        "SDM", \
        "TSDM", \
        "AT", \
        "NOP", \
        "AUC", \
        "DT", \
        "RT", \
        ])


def folder_retriever(folder):                               # input relative path "folder"
    # returns the "file" names as a list within the folder

    # makes a list of the files within "folder"
    folders_in_master = os.listdir(folder)
    folder_name = []                                # makes an empty list to store the folder name

    for i in folders_in_master:                     # iterates through each file within the "folder"
        # makes a string variable of the path to locate the "folder" and "files" within
        path_maker = "%s/%s" % (folder, i)
        if ("._" not in path_maker) and ("IONO" not in path_maker) and ("Iono" not in path_maker) and ("iono" not in path_maker):                   # excludes hidden files from folder
        # if "._"  not in path_maker:                  # excludes hidden files from folder
            folder_name.append(path_maker)             # fills the empty list with "files"

    return folder_name


# test code to see if function works
files_seperator = folder_retriever(input_main_folder)
# print(files_seperator)

z = 0
for folder_name in files_seperator:
    if input_folder_name == folder_name:
        # print(input_folder_name)                                                        # prints the extracted folder name
        input_folder_name = folder_name
        break
    z += 1



file_location_relative_to_main_folder = folder_retriever(files_seperator[z])             # prints the path of the "file" within the "folder"
# print(file_location_relative_to_main_folder)
print(len(file_location_relative_to_main_folder))

# exit()
y = 0
graphvisualiser = 0

# for z in file_location_relative_to_main_folder:

#11111999

for y in file_location_relative_to_main_folder:                               # this for loop finds the file (input_trial_run) in the folder (input_growth_factor)   

    cell_type = None                                                            # Initialize growth factor used on cell (GDNF vs NGF)
    sex = None                                                                  # Initialize sex of the animal used (Male or Female)
    serum_type = None                                                           # Initialize serum_type

    file_1 = str(y)
    newfile= file_1[(file_1.rfind("/")+1):]

    if "GDNF" in str(newfile):                                              #labels the columns to what the cell types are
        cell_type = "GDNF"

    if "NGF" in str(newfile):                                               #labels the columns to what the cell types are
        cell_type = "GDNF"
    


    start_num = 1                                                               # this block codes for the condition (COND) of the file analysed
    end_num = 20  # You can adjust this range as needed

    # Initialize serum_type


    # Iterate over the range of numbers
    for num in range(start_num, end_num + 1):
        # Construct the condition string
        condition_str = f"COND{num}"
        
        # Check if the condition string is present in newfile
        if condition_str in str(newfile):
            # If present, assign the condition string to serum_type
            serum_type = condition_str
            # Break out of the loop since we found the condition
            break

    # Check if serum_type was assigned
    if serum_type is not None:
        print(f"Condition found: {serum_type}")
    else:
        print("No condition found.")


    if "M" in str(newfile):                                              #labels the columns to what the cell types are
        sex = "M"
    else:
        sex = "F"



    # if "COND1" in str(newfile):                                             #labels the columns to what serum was added (the condition)
    #     serum_type = "COND1"

    # if "COND1" in str(newfile):                                             #labels the columns to what serum was added (the condition)
    #     serum_type = "COND1"



    graphvisualiser += 1

    column_importer = pd.read_csv(y)
    # print(column_importer)
    # exit()


    column_name = column_importer.columns.tolist()                                               # extracts the columns from the excel document
    intensity_mean = []                                                                          # makes a new empty list


    for i in column_name:
        # allows certain channels to be removed
        if (("Channel1" in i) and ("IntensityMeanThrs!!R" in i)) or "Relative Time!!R" in i:
        # if (("Channel1_R1" in i) and ("IntensityMeanThrs!!R" in i)) or "Relative Time!!R" in i:       # allows certain ROI to be graphed to be removed
        # if (("Channel2" in i) and ("IntensityMeanThrs!!R" in i)) or "Relative Time!!R" in i:          # allows certain channels to be removed
        # if (("Channel3" in i) and ("IntensityMeanThrs!!R" in i)) or "Relative Time!!R" in i:          # allows certain channels to be removed
        # if "IntensityMeanThrs!!R" in i or "Relative Time!!R" in i:
            intensity_mean.append(i)
    # print(column_importer[intensity_mean])                        # extracts the time column and all columns with the word "IntensityMeanThrs!!R"

    time = []                                                       # creates a time variable on floats and removes the s from the first row
    original_time_column = column_importer[intensity_mean[0]]
    for i in original_time_column[1:]:
        time.append(float(i))
        
    number_of_curves = 0
    non_responder = 0

    
    for j in column_importer[intensity_mean[1:-1]]:                                   # drops the first column as this is time, drops the last column as that is the control
    
        number_of_curves += 1 
        # print(j)
        
        intensity = []
        original_intensity_column = column_importer[j]                                # creates a intensity variable on floats and removes the s from the first row
        for i in original_intensity_column[1:]:
            intensity.append(float(i))
        
        
        gaussianoriginal = gaussian_filter1d(intensity, 0.0000001)                    # original data zeroed (no smoothing)
        # plt.plot(time, gaussianoriginal - gaussianoriginal[0])                      # plots the original data with translation to 0
        # plt.plot(time, gaussianoriginal)                                            # plots the original data (raw)


        gaussian_smoothed_curved = gaussian_filter1d(intensity, 6)                    # gaussian curve smoothing function
        zeroed_curves = gaussian_smoothed_curved - gaussian_smoothed_curved[0]        # zeros the curves (with smoothing)
        first_derivative = np.gradient(zeroed_curves)                                 # gets the first derivative of the zeroed_curves
        second_derivative = np.gradient(first_derivative)                             # gets the second derivative of the zeroed_curves
        

        maximum_fluorescence_intensity_position = np.squeeze(np.argwhere(zeroed_curves == max(zeroed_curves)))                                 # finds the time (nearest value) of the peak intensity
        maximum_fluorescence_intensity = zeroed_curves[maximum_fluorescence_intensity_position]                                                                     # finds the maximum point of the graph (zeroed_curves)
        time_to_maximum_fluorescence_intensity = time[maximum_fluorescence_intensity_position]                                                                           # finds the time the maximum point occured at of the graph (zeroed_curves)
        
        
        first_derivative_maximum_position = np.squeeze(np.argwhere(first_derivative == max(first_derivative)))               # finds the datapoint of the maximal velocity
        first_derivative_maximum = first_derivative[first_derivative_maximum_position]                                                   # finds the maximum velocity of the graph (zeroed_curves)
        time_to_first_derivative_maximum = time[first_derivative_maximum_position]                                                          # find the time the max velocity occurs at (maximal slope)
        
        
        second_derivative_maximum_position = np.squeeze(np.argwhere(second_derivative == max(second_derivative)))         # finds the datapoint of the maximal acceleration
        second_derivative_maximum = second_derivative[second_derivative_maximum_position]                                          # finds the maximum acceleraiton of the graph (zeroed_curves)
        time_to_second_derivative_maximum = time[second_derivative_maximum_position]                                                  # finds the time of the max acceleration (activaiton point)
        second_derivative_max_checker = find_peaks(second_derivative, prominence = 0.1)
        # print(second_derivative_max_checker)
        # print(find_peaks(zeroed_curves))
        location_of_acceleration_peaks = second_derivative_max_checker[0]
        # print(location_of_acceleration_peaks)
        # exit()

        for x in zeroed_curves[maximum_fluorescence_intensity_position::]:                                          
            if x < maximum_fluorescence_intensity/2:
                break
        decay_to_half_max_position = np.squeeze(np.argwhere(zeroed_curves == x))                                                # finds the position of the 50% decay point after the maximum
        decay_to_half_max_intensity = zeroed_curves[decay_to_half_max_position]
        decay_to_half_max_time = time[decay_to_half_max_position]


        area_under_curve = trapezoid(np.maximum(zeroed_curves, 0.0), time)                                                # calculates the area under the curve after translation and not accounting for values below the X-axis


        for serum_add in time:                                                                                   # creates a point for when the serum was added
            if serum_add > 10:
                break
        
        
        decay_time = decay_to_half_max_time - time_to_maximum_fluorescence_intensity                            #this needs to go before the next line as the time_to_maximum variable is changed
        time_to_maximum_fluorescence_intensity = time_to_maximum_fluorescence_intensity - serum_add
        first_derivative_maximum = first_derivative_maximum/float(time[1])
        time_to_first_derivative_maximum = time_to_first_derivative_maximum - serum_add
        second_derivative_maximum = second_derivative_maximum/float(time[1]) 
        time_to_second_derivative_maximum = time_to_second_derivative_maximum - serum_add
        activation_time = time_to_maximum_fluorescence_intensity - time_to_second_derivative_maximum
        

        peaks_identifier = find_peaks(zeroed_curves, height = zeroed_curves[second_derivative_maximum_position],
                                    distance = 5, prominence = 7)  # this gets the number of peaks
        location_of_peaks = peaks_identifier[0]                                                     # extracts the location (in the list) of the peaks
        # print(location_of_peaks)
        number_of_peaks = len(peaks_identifier[0])                                                  # prints the number of peaks
                
        # for i in location_of_peaks:
        #     time[i]                                                                                 # marks the time of the peak
        #     # print("time of peaks: " + str(time[i]))                                               # prints the time of the peaks
        #     zeroed_curves[i]                                                                        # marks the intensity of the peak   
        #     # print("intensity of peaks: " + str(zeroed_curves[i]))                                 # prints the intensity of the peaks

        plt.plot(time, zeroed_curves, color = 'green')

        run_time = time[-1]

        # if (time_to_maximum_fluorescence_intensity <= 0) or (time_to_first_derivative_maximum <= 0) or (time_to_second_derivative_maximum <= 0) \
        #     or (activation_time <= 0) or (maximum_fluorescence_intensity < 7) or (second_derivative_maximum < 0.05) or (number_of_peaks <= 0) or (decay_to_half_max_intensity > (maximum_fluorescence_intensity/2)):        # this is exclusion criteria for the curves if they are false
            
        #     # zeroed_curves = [0 for y in range(len(time))]
        #     # first_derivative = [0 for y in range(len(time))]
        #     # second_derivative = [0 for y in range(len(time))]

        #     maximum_fluorescence_intensity = None                                                                
        #     time_to_maximum_fluorescence_intensity = None   
        #     first_derivative_maximum = None                                            
        #     time_to_first_derivative_maximum = None
        #     second_derivative_maximum = None
        #     time_to_second_derivative_maximum = None
        #     number_of_peaks = None
        #     activation_time = None
        #     area_under_curve = None
        #     decay_time = None
        #     run_time = None
        #     non_responder += 1


        #     # plt.plot(time, first_derivative, label = "f'(x)")                                                      # graphs the first derivative of the graph (zeroed_curves)
        #     # plt.plot(time, second_derivative, color = 'blue')                                                 


        #     plt.plot(time, zeroed_curves, color = 'red') 
        # else:                                                                 ###the line above and below this can be turned on to show what has been included and exclused
        #     plt.plot(time, zeroed_curves, color = 'green')




        # print(
        #     "ROI: " + j,\
        #     "\nMaximum Fluorescence Intensity: " + str(maximum_fluorescence_intensity), \
        #     "\nTime to Maximum Fluorescence Intensity: " + str(time_to_maximum_fluorescence_intensity), \
        #     "\nFirst Derivative Maximum: " + str(first_derivative_maximum), \
        #     "\nTime to First Derivative Maximum: " + str(time_to_first_derivative_maximum), \
        #     "\nSecond Derivative Maximum: " + str(second_derivative_maximum), \
        #     "\nTime to Second Derivative Maximum: " + str(time_to_second_derivative_maximum),\
        #     "\nActivation Time: " + str(activation_time), \
        #     "\nNumber of peaks: " + str(number_of_peaks), \
        #     "\nArea Under Curve: " + str(area_under_curve), \
        #     "\n"
        #     )


        plt.subplot(round((len(file_location_relative_to_main_folder) + 1)/3), 3, graphvisualiser)
        plt.title(y)
        plt.xlabel("Time (Seconds)")
        plt.ylabel("Intensity")
        # plt.xticks(np.arange(0, 200, 50))
        # plt.legend()
        # plt.grid()

        
        
        # plt.plot(time, zeroed_curves, label = "f(x)")                                                                            # graphs the zeroed curves   
        # plt.plot(time, first_derivative, label = "f'(x)")                                                      # graphs the first derivative of the graph (zeroed_curves)
        # plt.plot(time, second_derivative, label = "f''(x)")                                                    # graphs the second derivative of the graph (zeroed_curves)

        with open(str(input_folder_name +  '_feature_extraction.csv'), 'a', newline='') as f:              # appends the made file with the data extraction
            thewriter = csv.writer(f)
            thewriter.writerow([
                y, \
                cell_type, \
                serum_type, \
                sex, \
                number_of_curves, \
                maximum_fluorescence_intensity, \
                time_to_maximum_fluorescence_intensity, \
                first_derivative_maximum, \
                time_to_first_derivative_maximum, \
                second_derivative_maximum, \
                time_to_second_derivative_maximum, \
                activation_time, \
                number_of_peaks, \
                area_under_curve, \
                decay_time, \
                run_time, \
                ])


    # print ("number of curves" + str(number_of_curves))
    # print ("number of non-responders" + str(non_responder))
    percentage_of_responders = float(((number_of_curves - non_responder)/number_of_curves)*100)
    # print ("Percentage of responders: " + str(percentage_of_responders))
    

plt.show()
