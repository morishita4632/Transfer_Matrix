import glob, os
lattice = "Xsquare"
Js = "2_2_1_1"

path = "out/" + lattice + '/' + Js + '/'
wfile = path + Js + "_all.txt"
if os.path.isfile(wfile): print("already merged")

file_list = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

with open(wfile, mode='w') as wf:
    for file in file_list:
        if file[-4:] == ".txt": continue
        rfile = path + file
        with open(rfile) as rf:
            line = 0
            for s_line in rf:
                line += 1
                if line % 9 == 3:
                    wf.write(s_line)