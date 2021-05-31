import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--overview", help = "Input result file")
parser.add_argument("-t", "--time", help = "Output overview file")
args, unknown = parser.parse_known_args()

def readjson(data_path):
    with open(data_path) as f:
        content = f.readlines()

    data = {}
    data["Step_3"] = content[0][:-2]
    data["Precision"] = float(content[1][:-2])
    data["Recall"] = float(content[2][:-2])
    data["Fmeasure"] = float(content[3][:-2])

    return data

data = readjson(args.overview)

with open(args.time, 'r+') as outfile:
    file = json.load(outfile)
    data["Execution"] = float(file["Step_1_1"]) + float(file["Step_1_2"]) + float(file["Step_1_3"]) + float(file["Step_2_1"]) + float(file["Step_2_2"]) + float(data["Step_3"])
    file.update(data)
    outfile.seek(0)
    json.dump(file, outfile)

