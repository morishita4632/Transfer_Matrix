import glob
import os
lattice = "Xsquare"
rfile = "out/Xsquare/slurm-106078.out"

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
        if line % 8 == 3:
            with open(dir + '/' + Js_str + "_even_all.txt", mode='a') as wf:
                wf.write(s_line)
        if line % 8 == 7:
            with open(dir + '/' + Js_str + "_odd_all.txt", mode='a') as wf:
                wf.write(s_line)
