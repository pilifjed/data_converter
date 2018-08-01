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

def float_convertable(value):
    """Checks if value is converable to float"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def parse_calibration_file(file_name):
    """Reads text digi file and stores it in a map structure"""
    CHANELS_NO = 16
    SAMPIC_NO = 2
    VALUES_NO = CHANELS_NO*SAMPIC_NO
    CELLS_NO = 64

    def parse_line(line):
        return list(map(float,line.split(" ")[:-1]))
    def end_of_channel(line_counter):
        return (line_counter % CELLS_NO == 0)

    with io_open(file_name,"r",errors="replace") as text_file:
        first = True
        line_counter = 1
        temp_cells = []
        output = {"values":[]}
        temp_val = {"db": -1, "sampic": -1, "chanel": -1, "cells": []}
        for line in text_file:
            if(float_convertable(line.split(" ")[0])):
                if(first):
                    first = False
                    output["voltage"] = parse_line(line)
                else:
                    temp_cells.append(parse_line(line))
                    if(end_of_channel(line_counter)):
                        try:
                            db = re.search('_(\d+)',file_name)
                            temp_val["db"] = int(db.group(1))
                        except(AttributeError):
                            print("DB not specified in file name")
                        temp_val["sampic"] = 1 - ((line_counter - 1) // CELLS_NO) // CHANELS_NO
                        temp_val["chanel"] = ((line_counter - 1) // CELLS_NO) % 16
                        temp_val["cells"] = temp_cells
                        output["values"].append(temp_val.copy())
                        temp_cells = []
                    line_counter+=1

        if(VALUES_NO != len(output["values"])):
            raise FileParseException('Corrupt file or wrong CHANELS_NO, SAMPIC_NO or CELLS_NO values specified. Expected CHANELS_NO*SAMPIC_NO = '
                                     + str(VALUES_NO) +', got '+ str(len(output["values"])))
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

def handle_optional_parameters(data, file_name, output_file_name, fit_from, fit_to):
    """Sets default values of optional parameters if unset"""
    if((fit_from is None) or (float(fit_from) < 0)):
        fit_from=0

    max_fit_to = max(data["voltage"])
    if((fit_to is None) or (float(fit_to) > max_fit_to)):
        fit_to=max_fit_to

    if(output_file_name is None):
        output_file_name = file_name[:file_name.rfind('.')] + '.cal.json' # If no other filename specified save in the same location as input file
    else:
        if output_file_name[-1] != "/":
            output_file_name += "/"
        output_file_name += file_name[file_name.rfind('/')+1:file_name.rfind('.')] + '.cal.json'
    return output_file_name, float(fit_from), float(fit_to)

def draw_plots(x, all_y, fit, chanel_to_plot):
    """Draws plots"""
    plt.figure(1)
    plt.title("Chanel "+ chanel_to_plot + " cells")
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


def convert_calibration_file(file_name, function, output_file_name=None, fit_from=None, fit_to=None, chanel_to_plot=None):
    """Main calibration function"""
    data = parse_calibration_file(file_name)
    output_file_name, fit_from, fit_to = handle_optional_parameters(data, file_name, output_file_name, fit_from, fit_to)

    fit = ROOT.TF1("fit",function, fit_from, fit_to)

    if(chanel_to_plot is not None):
        chanel_no = chanel_to_plot
        draw_plots(data["voltage"], data["values"][int(chanel_no)]["cells"], fit, chanel_no)

    all_values = data["values"]
    for value in all_values:
        conv_value = []
        for cell in value["cells"]:
            conv_value.append(process_cell(data["voltage"], cell, fit))
        value["cells"] = conv_value
    del data["voltage"]
    data["parameters"] = data.pop("values")
    data["formula"] = function

    with open(output_file_name, 'w') as outfile:
        json.dump(data, outfile)

def main():
    parser = argparse.ArgumentParser(description='Converts text calibration data to .cal.JSON format,by fitting function to the given data and storing its parameters.')
    parser.add_argument("input_file",action="store",help="set input file name")
    parser.add_argument("function",action="store",help="set function to fit i.e. \'[0]*x+[1]\'")
    parser.add_argument("-o", "--output", action="store", metavar= "PATH/TO/FILE", help="set path where file schould be stored")
    parser.add_argument("-f", "--fit_from", action="store", metavar= "LEFT_BOUND", help="set left bound of fit values interval")
    parser.add_argument("-t", "--fit_to", action="store", metavar= "RIGHT_BOUND", help="set right bound of fit values interval")
    parser.add_argument("-p", "--plot", action="store", metavar= "CHANEL_NO", help="enable plotting mode, takes number of channel to plot")
    args = parser.parse_args()
    convert_calibration_file(args.input_file, args.function, args.output, args.fit_from, args.fit_to, args.plot)

if __name__ == '__main__':
    main()
