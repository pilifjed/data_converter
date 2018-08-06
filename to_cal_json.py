import json
from io import open as io_open
import argparse
import re
import matplotlib.pyplot as plt
import numpy as np

import ROOT
from array import array


class FileParseException(Exception):
    pass

class FileMergeException(Exception):
    pass

def float_convertable(value):
    """Checks if value is converable to float"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def parse_calibration_file(file_name):
    """Reads text digi file and stores it in a map structure"""
    CHANNELS_NO = 16
    SAMPIC_NO = 2
    VALUES_NO = CHANNELS_NO*SAMPIC_NO
    CELLS_NO = 64

    def parse_line(line):
        return list(map(float,line.split(" ")[:-1]))
    def end_of_channel(line_counter):
        return (line_counter % CELLS_NO == 0)

    with io_open(file_name,"r",errors="replace") as text_file:
        first = True
        line_counter = 1
        temp_cells = []
        output = {"values":{}}
        temp_val = {}

        try:
            db = re.search('_db(\d+)',file_name, re.IGNORECASE)
            temp_val["db"] = int(db.group(1))
        except AttributeError:
            raise FileParseException("DB not specified in file name")

        output["values"][temp_val["db"]] = []
        for line in text_file:
            if(float_convertable(line.split(" ")[0])):
                if(first):
                    first = False
                    output["voltage"] = parse_line(line)
                else:
                    temp_cells.append(parse_line(line))
                    if(end_of_channel(line_counter)):
                        temp_val["time_offset"] = 0
                        temp_val["sampic"] = 1 - ((line_counter - 1) // CELLS_NO) // CHANNELS_NO
                        temp_val["channel"] = ((line_counter - 1) // CELLS_NO) % 16
                        temp_val["cells"] = temp_cells
                        output["values"][temp_val["db"]].append(temp_val.copy())
                        temp_cells = []
                    line_counter+=1

        if(VALUES_NO != len(output["values"][next(iter(output["values"]))])):
            raise FileParseException('Corrupt file or wrong CHANNELS_NO, SAMPIC_NO or CELLS_NO values specified. Expected CHANNELS_NO*SAMPIC_NO = '
                                     + str(VALUES_NO) +', got '+ str(len(output["values"][next(iter(output["values"]))])))
    return output

def process_cell(voltage, cell, fit):
    """Fits given TF1 to the cell, and stores its parameters"""
    x = array('d',sorted(voltage))
    y = array('d',sorted(cell))
    graph = ROOT.TGraph(len(x),x,y)
    parameters = graph.Fit(fit, 'SQR')
    output = []
    for i in range(0,fit.GetNpar()):
        output.append(parameters.Value(i))
    return output

def handle_optional_parameters(parsed_data, fit_from, fit_to):
    """Sets default values of optional parameters if unset"""
    if((fit_from is None)): # """or (float(fit_from) < 0)"""
        fit_from=0
    if((fit_to is None)): # """or (float(fit_to) > max_fit_to)"""
        max_fit_to = min(max((float(x) for x in data["voltage"])) for file_name, data in parsed_data)
        fit_to=max_fit_to
    return  float(fit_from), float(fit_to)



def draw_plots(x, all_y, fit, channel_to_plot):
    """Draws plots"""
    plt.figure(1)
    plt.title("Channel "+ channel_to_plot + " cells")
    plt.xlabel("Voltage")
    plt.ylabel("Digital value")
    for y in all_y:
        plt.plot(x,y)
    process_cell(x, all_y[0], fit)
    fit_val = [fit.Eval(i) for i in x]
    residual_val = [abs(i-j) for i,j in zip(fit_val,all_y[0])]
    fig = plt.figure(2)
    plt.title("Values and Fit")
    plt.xlabel("Voltage")
    plt.ylabel("Digital value")
    plt.plot(x,all_y[0])
    plt.plot(x,fit_val)
    plt.legend(["Digital values","Fitted function"],bbox_to_anchor=(1, 1), loc=4, borderaxespad=0.)
    plt.figure(3)
    plt.title("Fit residual")
    plt.xlabel("Voltage")
    plt.ylabel("Residual")
    plt.scatter(x,residual_val)
    plt.show()


def convert_calibration_file(data, function, fit_from=None, fit_to=None, channel_to_plot=None):
    """Main calibration function"""

    fit = ROOT.TF1("fit",function, fit_from, fit_to)

    db_no = next(iter(data["values"]))

    if(channel_to_plot is not None):
        channel_no = channel_to_plot
        draw_plots(data["voltage"], data["values"][db_no][int(channel_no)]["cells"], fit, channel_no)

    all_values = data["values"][db_no]
    for value in all_values:
        conv_value = []
        for cell in value["cells"]:
            conv_value.append(process_cell(data["voltage"], cell, fit))
        value["cells"] = conv_value
    del data["voltage"]
    data["parameters"] = data.pop("values")
    data["formula"] = function
    return data

def read_cal_json(file_name):
    with open(file_name) as file:
        data = json.load(file)
    return data

def merge(to_merge):
    merged_data = to_merge.pop()
    merged_keys = list(merged_data["parameters"].keys())
    for data in to_merge:
        for key in data["parameters"]:
            if key not in merged_keys:
                merged_keys.append(key)
                merged_data["parameters"][key] = data["parameters"][key]
            else:
                raise FileMergeException('Trying to merge files with the same digitizer board')
    return merged_data

def has_extension(extension, file_name):
    return re.search(extension, file_name, re.IGNORECASE) is not None

def main():
    parser = argparse.ArgumentParser(description='Converts text calibration data to .cal.json format,by fitting function to the given data and storing its parameters.')
    parser.add_argument("input_files", action="store", nargs="+", help="set input file names")
    parser.add_argument("-f", "--function", default="pol1",action="store",help="set function to fit i.e. \'[0]*x+[1]\'")
    parser.add_argument("-l", "--fit_from", action="store", metavar="LEFT_FIT_BOUND", help="set left bound of fit values interval")
    parser.add_argument("-r", "--fit_to", action="store", metavar="RIGHT_FIT_BOUND", help="set right bound of fit values interval")
    parser.add_argument("-p", "--plot",action="store", metavar="CHANNEL_TO_PLOT_NO", help="enable plotting during conversion")
    parser.add_argument("-m", "--merge", action="store", metavar="OUTPUT_FILE_NAME", help="enable merge mode")
    args = parser.parse_args()

    parsed_data = []
    converted_data = []
    output_data = []



    for file_name in args.input_files:
        if has_extension(".cal.json", file_name):
            data = read_cal_json(file_name)
            function = data["formula"]
            if args.merge and function != args.function:
                print(file_name + "didn't match function specified for conversion, file won't be included in merged output file")
            else:
                output_data.append((file_name, data))
        else:
            try:
                data = parse_calibration_file(file_name)
                parsed_data.append((file_name, data))
            except Exception as e:
                print("File " + file_name + " skipped. An error ocurred during parse process: " + str(e))

    if parsed_data != []:
        fit_from, fit_to = handle_optional_parameters(parsed_data, args.fit_from, args.fit_to)

    for file_name, data in parsed_data:
        try:
            new_data = convert_calibration_file(data, args.function, fit_from, fit_to, args.plot)
            converted_data.append((file_name, new_data))
            args.plot = None # only first file will be plotted
        except Exception as e:
            print("File " + file_name + " skipped. An error ocurred during conversion process: " + str(e))

    for file_name, data in converted_data:
        output_file_name = file_name[:file_name.rfind(".")] + ".cal.json"
        output_data.append((output_file_name, data))

    if args.merge is not None:
        merge_file_name = args.merge
        if not has_extension(".cal.json", args.merge):
            merge_file_name += ".cal.json"
        merged_data = merge([data for file_name, data in output_data])
        output_data = [(merge_file_name, merged_data)]

    for output_file_name, data in output_data:
        with open(output_file_name, 'w') as outfile:
            json.dump(data, outfile)

if __name__ == '__main__':
    main()
