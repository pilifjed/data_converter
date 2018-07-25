import json
from io import open as io_open
import argparse

import ROOT
from array import array


class FileParseException(Exception):
    pass

def float_convertable(value):
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
    parameters = graph.Fit(fit, 'S')
    output = []
    for i in range(0,2):
        output.append(parameters.Value(i))
    return output

def handle_optional_parameters(data, file_name, output_file_name, fit_from, fit_to):
    """Sets default values of optional parameters if unset"""
    if(fit_from is None or fit_from < 0):
        fit_from=0

    max_fit_to = len(data["voltage"])
    if(fit_to is None or fit_to > max_fit_to):
        fit_to=max_fit_to

    if(output_file_name is None):
        output_file_name = file_name[:file_name.rfind('.')] + '.cal.json' # If no other filename specified save in the same location as input file
    else:
        output_file_name += '.cal.json'
    return output_file_name, fit_from, fit_to


def convert_calibration_file(file_name, function, output_file_name=None, fit_from=None, fit_to=None):
    """ """
    data = parse_calibration_file(file_name)
    output_file_name, fit_from, fit_to = handle_optional_parameters(data, file_name, output_file_name, fit_from, fit_to)

    fit = ROOT.TF1("fit",function, fit_from, fit_to)
    all_values = data["values"]
    for value in all_values:
        conv_value = []
        for cell in value["cells"]:
            conv_value.append(process_cell(data["voltage"], cell, fit))
        value["cells"] = conv_value

    del data["voltage"]
    data["function"] = function

    with open(output_file_name, 'w') as outfile:
        json.dump(data, outfile)

def main():
    parser = argparse.ArgumentParser(description='Converts text calibration data to .cal.JSON format,by fitting function to the given data and storing its parameters.')
    parser.add_argument("input_file",action="store",help="Input file name")
    parser.add_argument("function",action="store",help="Function to fit i.e. \"[0]*x+[1]\"")
    parser.add_argument("-o", "--output", action="store", help="Path to the location where output will be stored")
    parser.add_argument("-f", "--fit_from", action="store", help="Left bound of fit values interval")
    parser.add_argument("-t", "--fit_to", action="store", help="Right bound of fit values interval")
    args = parser.parse_args()
    convert_calibration_file(args.input_file, args.function, args.output, args.fit_from, args.fit_to)

if __name__ == '__main__':
    main()
