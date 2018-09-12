# Usage: python eval_dataset_eye_tracking.py <path to text file to be analysed> [<location to store analysed file>]
import sys, os, re, importlib, statistics, math, operator, csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, integrate, interpolate
from numpy import NaN, Inf, arange, isscalar, asarray, array

def average_x(x):
    try:
        return str(statistics.mean(list(map(float, x))))
    except:
        # print("exception in average")
        return str(0)

def sum_x(x):
    try:
        return str(sum(list(map(float, x))))
    except:
        # print("exception in sum")
        return str(0)

def first_x(x):
    try:
        return x[0]
    except:
        # print("exception in first")
        return "."

def count_x(x):
    try:
        return str(len(x))
    except:
        # print("exception in len")
        return str(0)

def float_x(x):
    y = []
    for i in range(len(x)):
        try:
            y.append(float(x[i]))
        except:
            pass
            # y.append(0)
    return y

headers = [('EYE_TRACKED', first_x), ('RESOLUTION_X', average_x), ('RESOLUTION_Y', average_x)]
eye_tracked = [('_ACCELERATION_X', average_x), ('_ACCELERATION_Y', average_x), ('_FIX_INDEX', average_x), ('_GAZE_X', average_x), ('_GAZE_Y', average_x), ('_INTEREST_AREAS', first_x), ('_INTEREST_AREA_DATA', first_x), ('_INTEREST_AREA_DISTANCE', average_x), ('_INTEREST_AREA_ID', first_x), ('_INTEREST_AREA_LABEL', first_x), ('_INTEREST_AREA_PIXEL_AREA', average_x), ('_IN_BLINK', average_x), ('_IN_SACCADE', average_x), ('_PUPIL_SIZE', average_x), ('_SACCADE_INDEX', average_x), ('_VELOCITY_X', average_x), ('_VELOCITY_Y', average_x)]
details = [('SAMPLE_BUTTON', average_x), ('SAMPLE_INDEX', average_x), ('SAMPLE_INPUT', first_x), ('SAMPLE_MESSAGE', first_x), ('TARGET_ACCELERATION_X', average_x), ('TARGET_ACCELERATION_Y', average_x), ('TARGET_VELOCITY_X', average_x), ('TARGET_VELOCITY_Y', average_x), ('TIMESTAMP', average_x), ('TRIAL_LABEL', first_x), ('TRIAL_START_TIME', first_x), ('VARIABLE_KEY', first_x), ('VARIABLE_TIME', average_x), ('height', average_x)]

# print("command line arguments:", sys.argv[:])
program_name = sys.argv[0]
filepath = sys.argv[1]
storage_location = ""
plot_graph = False
if len(sys.argv) > 2:
    storage_location = sys.argv[2]
if not re.search("/$", storage_location):
    storage_location = storage_location + "/"
# print("program_name:", program_name)
# print("filepath:", filepath)
if os.path.exists(filepath):
    filename  = os.path.splitext(os.path.basename(filepath))[0] # filename without extension
    extension = os.path.splitext(os.path.basename(filepath))[1] # extension of the file along with separater
    csvfile = open(filepath, "r")
    # print()
    # print(filename)

    if extension != ".csv":
        print("not a .csv file...")
        sys.exit()

    rows = list(csv.reader(csvfile))
    # print(len(rows[2]))
    # sys.exit()
    first_line = rows[0]
    second_line = rows[1]
    # print("first_line:", first_line)
    # print("second_line:", second_line)
    default_split = True
    sampling_rate = 500
    column_count = 0
    drawing = False
    if first_line[0] == "question":
        first_line = first_line[1:]
        second_line = second_line[1:]
        # print(second_line)
        second_line = [ int(math.floor(float(num)*sampling_rate)) for num in second_line ]
        default_split = False
        rows = rows[2:]
        column_count = len(rows[2])
        if first_line[-1] == "drawing":
            drawing = True
            del first_line[-1]
        # print(os.path.basename(filename))
        # print(first_line)
        # print(second_line)
        new_second_line = []
        for que in first_line:
            # print(que)
            search_string = "p" + que + " start"
            # print(search_string)
            start_index = [i for i, row in enumerate(rows) if search_string in row]
            # print(start_index)
            try:
                new_second_line.append(start_index[0])
            except:
                new_second_line = second_line
                break;
            search_string = "p" + que + " end"
            # print(search_string)
            end_index = [i for i, row in enumerate(rows) if search_string in row]
            # print(end_index)
            try:
                new_second_line.append(end_index[0])
            except:
                new_second_line = second_line
                break;
        # print(new_second_line)
        second_line = new_second_line
        # print(second_line)
        # sys.exit()
    else:
        column_count = len(rows[0])
    # print(len(rows[0]))
    # print(rows[0])
    eye_index = rows[0].index('EYE_TRACKED')
    # print(eye_index)
    eye = rows[1][eye_index]
    # print(eye)
    if eye == 'Right':
        eye_tracked = list(map(lambda x: ("RIGHT" + x[0], x[1]), eye_tracked))
    else:
        eye_tracked = list(map(lambda x: ("LEFT" + x[0], x[1]), eye_tracked))
    # print(all(map(lambda x: len(x) == len(rows[0]), rows[1:])))
    # print([all(list(rows[0].index(label[0]) for label in headers))], [all(list(rows[0].index(label[0]) for label in eye_tracked))], [all(list(rows[0].index(label[0]) for label in details))])

    fig = plt.gcf()
    fig.set_size_inches(25, 9, forward=True)

    interest_area_index = rows[0].index(eye.upper() + "_INTEREST_AREAS")
    start_pointer, end_pointer = 1, 1

    analysed_data = []
    temp_analysed_data = []

    for label in headers:
        temp_analysed_data.append(label[0])
    for label in eye_tracked:
        temp_analysed_data.append(label[0])
    for label in details:
        temp_analysed_data.append(label[0])
    analysed_data.append(temp_analysed_data)

    if default_split:
        while start_pointer < len(rows):
            while end_pointer < len(rows) and rows[start_pointer][interest_area_index] == rows[end_pointer][interest_area_index]:
                end_pointer = end_pointer + 1
            filtered_data = []
            for label in headers:
                try:
                    # print(label)
                    ind = rows[0].index(label[0])
                    # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                    filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                except:
                    filtered_data.append("ERROR: NA")
            for label in eye_tracked:
                try:
                    # print(label)
                    ind = rows[0].index(label[0])
                    # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                    filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                except:
                    filtered_data.append("ERROR: NA")
            for label in details:
                try:
                    # print(label)
                    ind = rows[0].index(label[0])
                    # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                    filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                except:
                    filtered_data.append("ERROR: NA")
            analysed_data.append(filtered_data)
            start_pointer = end_pointer
    else:
        current_index = 0
        summary_eye_tracking = {"total_pupil_size" : [], "pupil_size_samples" : [], "details" : ""}
        summary_eye_tracking["details"] = summary_eye_tracking["details"] + "\n" + str(filename) + "\n"
        while current_index < len(first_line):
            # Question wise analysis
            # print("##############")
            # print(current_index)
            filtered_data = []
            filtered_data = ["question", first_line[current_index]]
            for i in range(column_count - 2):
                filtered_data.append(" ")
            analysed_data.append(filtered_data)
            summary_eye_tracking["details"] = summary_eye_tracking["details"] + "Question: " + str(first_line[current_index]) + "\n"
            # print("Question:", first_line[current_index])
            try:
                ind = rows[0].index(eye.upper() + "_PUPIL_SIZE")
                # print(list(map(float, [rows[i][ind] for i in range(second_line[2*current_index], second_line[2*current_index + 1])])))
                summary_eye_tracking["total_pupil_size"].append(sum(float_x([rows[i][ind] for i in range(second_line[2*current_index], second_line[2*current_index + 1])])))
                summary_eye_tracking["pupil_size_samples"].append(len(float_x([rows[i][ind] for i in range(second_line[2*current_index], second_line[2*current_index + 1])])))
                summary_eye_tracking["details"] = summary_eye_tracking["details"] + "Average pupil size: " + str(round(summary_eye_tracking["total_pupil_size"][-1]/summary_eye_tracking["pupil_size_samples"][-1], 3)) + "\n"
                # print("Average pupil size:", round(summary_eye_tracking["total_pupil_size"][-1]/summary_eye_tracking["pupil_size_samples"][-1], 3))
            except:
                summary_eye_tracking["details"] = summary_eye_tracking["details"] + "Average pupil size: NA" + "\n"
                # print("Average pupil size:", round(summary_eye_tracking["total_pupil_size"][-1]/summary_eye_tracking["pupil_size_samples"][-1], 3))
            start_pointer = second_line[2*current_index]
            end_pointer   = second_line[2*current_index]
            flag = False
            eye_tracking_data = {"x": [], "interest_area":[], "pupil_size":[], "saccade_length":[]}
            x_counter = 0
            while start_pointer < len(rows) and start_pointer < second_line[2*current_index + 1]:
                # Same interest ID
                flag = True
                # print("++++++++++++++")
                # print(start_pointer)
                while end_pointer < len(rows) and end_pointer <= second_line[2*current_index + 1] and rows[start_pointer][interest_area_index] == rows[end_pointer][interest_area_index]:
                    end_pointer = end_pointer + 1
                # print(end_pointer)
                # print("++++++++++++++")
                filtered_data = []
                for label in headers:
                    try:
                        # print(label)
                        ind = rows[0].index(label[0])
                        # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                        filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                    except:
                        filtered_data.append("ERROR: NA")
                for label in eye_tracked:
                    try:
                        # print(label)
                        ind = rows[0].index(label[0])
                        # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                        filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                    except:
                        filtered_data.append("ERROR: NA")
                for label in details:
                    try:
                        # print(label)
                        ind = rows[0].index(label[0])
                        # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                        filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                    except:
                        filtered_data.append("ERROR: NA")
                analysed_data.append(filtered_data)
                # Creating array for graphs
                # if(re.match('.*\d+.*', rows[start_pointer][interest_area_index])):
                x_counter = x_counter + 1
                try:
                    timestamp_index = rows[0].index("TIMESTAMP")
                    trial_start_time_index = rows[0].index("TRIAL_START_TIME")
                    eye_tracking_data["x"].append(float(rows[start_pointer][timestamp_index]) - float(rows[start_pointer][trial_start_time_index]))
                except:
                    eye_tracking_data["x"].append(x_counter)
                eye_tracking_data["interest_area"].append(re.sub('[^0-9]', '' ,rows[start_pointer][interest_area_index]))
                # print(rows[start_pointer][interest_area_index])
                try:
                    ind = rows[0].index(eye.upper() + "_PUPIL_SIZE")
                    # print("*", ind)
                    # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                    eye_tracking_data["pupil_size"].append(float(average_x([rows[i][ind] for i in range(start_pointer, end_pointer)])))
                except:
                    eye_tracking_data["pupil_size"].append(-1)
                try:
                    ind = rows[0].index(eye.upper() + "_IN_SACCADE")
                    # print("*", ind)
                    # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                    eye_tracking_data["saccade_length"].append(float(average_x([rows[i][ind] for i in range(start_pointer, end_pointer)])))
                except:
                    eye_tracking_data["saccade_length"].append(-1)
                start_pointer = end_pointer
            # print(flag)
            # print("##############")
            # Plotting Graphs
            if(plot_graph):
                # print("*", len(eye_tracking_data["x"]))
                fig = plt.gcf()
                fig.set_size_inches(max(math.ceil(eye_tracking_data["x"][-1] - eye_tracking_data["x"][0])/(5*sampling_rate), 16), 9, forward=True)
                plt.plot(eye_tracking_data["x"], eye_tracking_data["pupil_size"])
                plt.scatter(eye_tracking_data["x"], eye_tracking_data["pupil_size"], c='r')
                plt.xlim(eye_tracking_data["x"][0], eye_tracking_data["x"][-1])
                plt.xticks(eye_tracking_data["x"], eye_tracking_data["interest_area"], fontsize=8)
                plt.xlabel("Interest Area", fontsize=8)
                plt.ylabel("Pupil Size", fontsize=8)
                plt.title("Question " + str(first_line[current_index]))
                plt.show()
                # plt.savefig(storage_location + filename + "_analysis_pupil_size(" + str(first_line[current_index]) + ").png", format="png", dpi=100)
                plt.gcf().clear()
                # sys.exit()
                # plt.plot(eye_tracking_data["x"], eye_tracking_data["saccade_length"])
                # plt.scatter(eye_tracking_data["x"], eye_tracking_data["saccade_length"], c='r')
                # plt.xticks(eye_tracking_data["x"], eye_tracking_data["interest_area"])
                # plt.xlabel("Interest Area")
                # plt.ylabel("Saccade Length")
                # plt.title("Question " + str(first_line[current_index]))
                # # plt.show()
                # plt.savefig(storage_location + filename + "_analysis_saccade_length(" + str(first_line[current_index]) + ").png", format="png", dpi=1000)
                # plt.gcf().clear()
                # sys.exit()
            current_index = current_index + 1
        if drawing:
            filtered_data = []
            filtered_data = ["drawing"]
            for i in range(column_count - 1):
                filtered_data.append(" ")
            analysed_data.append(filtered_data)
            while start_pointer < len(rows):
                while end_pointer < len(rows) and rows[start_pointer][interest_area_index] == rows[end_pointer][interest_area_index]:
                    end_pointer = end_pointer + 1
                filtered_data = []
                for label in headers:
                    try:
                        # print(label)
                        ind = rows[0].index(label[0])
                        # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                        filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                    except:
                        filtered_data.append("ERROR: NA")
                for label in eye_tracked:
                    try:
                        # print(label)
                        ind = rows[0].index(label[0])
                        # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                        filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                    except:
                        filtered_data.append("ERROR: NA")
                for label in details:
                    try:
                        # print(label)
                        ind = rows[0].index(label[0])
                        # print([rows[i][ind] for i in range(start_pointer, end_pointer)])
                        filtered_data.append(label[1]([rows[i][ind] for i in range(start_pointer, end_pointer)]))
                    except:
                        filtered_data.append("ERROR: NA")
                analysed_data.append(filtered_data)
                start_pointer = end_pointer
        # print()
        try:
            summary_eye_tracking["details"] = summary_eye_tracking["details"] + "\n" + "Average pupil size for entire crt: " + str(round(sum(summary_eye_tracking["total_pupil_size"])/sum(summary_eye_tracking["pupil_size_samples"]), 3))
            # print("Average pupil size for entire crt:", round(sum(summary_eye_tracking["total_pupil_size"])/sum(summary_eye_tracking["pupil_size_samples"]), 3))
        except:
            summary_eye_tracking["details"] = summary_eye_tracking["details"] + "\n" + "Average pupil size for entire crt: NA"
        print(summary_eye_tracking["details"])


    try:
        f = open(storage_location + filename + "_analysis" + extension, 'w+')
        writer = csv.writer(f)
        writer.writerows(analysed_data)
        f.close()
    except:
        print("Error in opening", storage_location + filename + "_analysis" + extension)
        print("Displaying on terminal")
        print('\n'.join(','.join(row) for row in analysed_data))

else:
    print("No such file exists")
