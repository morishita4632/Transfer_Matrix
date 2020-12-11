import glob
import os
lattice = "Xsquare"
suffix = "_odd"
rfile = "out/Xsquare/slurm.out"

line = 0
Js_str = ""
dir = ""

with open(rfile) as rf:
    for s_line in rf:
        line += 1
        if line % 52 == 1:
            Js_raw = [float(s) for s in s_line[5:-2].split(", ")]
            Js_str = '_'.join([str(int(J) if J.is_integer() else J)
                               for J in Js_raw])
            print(Js_str)
            dir = "out/" + lattice + '/' + Js_str
            if not os.path.isdir(dir):
                os.makedirs(dir)
        if line % 4 == 3:
            with open(dir + '/' + Js_str + "_all" + suffix + ".txt", mode='a') as wf:
                wf.write(s_line)
