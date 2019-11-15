#!/usr/bin/env python

### Usage: .\calcPseudoseq.py T0_CTL.txt T0.txt test.txt

## Project: Pseudo-Seq
## Purpose: Count differential peaks between CTL and CMC samples
## Created By: Gabriela Toomer
## Created Date: November 11, 2019
#Read the count data

from sys import argv, exit

try:
    xrange
except NameError:
    xrange = range

half_ws = 75 # window size is 150nt

def calc_sum_array(tlen, vec):
    sum_arr = [0.0] * tlen
    sum_arr[0] = vec[0]
    for i in xrange(1, tlen):
        sum_arr[i] = vec[i] + sum_arr[i - 1]
    return sum_arr
    
def determine_ws(pos, csum, tsum):
    [ws, wr_plus, wr_minus] = [0.0, 0.0, 0.0]

    # left half window
    ub = pos - 1
    lb = ub - half_ws
    if ub >= 0:
        if lb < 0:
            ws += ub + 1
            wr_plus += tsum[ub]
            wr_minus += csum[ub]
        else:
            ws += half_ws
            wr_plus += tsum[ub] - tsum[lb]
            wr_minus += csum[ub] - csum[lb]

    # right half window
    ub = min(pos + half_ws, tlen - 1)
    lb = pos
    ws += ub - lb
    wr_plus += tsum[ub] - tsum[lb]
    wr_minus += csum[ub] - csum[lb]

    ws = half_ws * 2 # the window size is fixed

    return [ws, wr_plus, wr_minus]
    

if len(argv) != 4:
    print("Usage: calcPseudoseq control.cnts treatment.cnts output")
    exit(-1)

with open(argv[1], "r") as fc:
    with open(argv[2], "r") as ft:
        with open(argv[3], "w") as fout:
            line_no = 0
            for line in fc:
                line_no += 1
                line2 = next(ft)

                # Parse data
                control = [float(x) for x in line.strip().split()]
                treatment = [float(x) for x in line2.strip().split()]
                tlen = len(control)
                assert tlen == len(treatment)
                # Normalize control
                sum_c = sum(control)
                if sum_c < 1e-8:
                    fout.write(" ".join(["0"] * tlen) + "\n")
                    continue
                factor = sum(treatment) / sum_c
                control = [x * factor for x in control]
                # Calculate psi scores
                ## Calculate sums 
                csum = calc_sum_array(tlen, control)
                tsum = calc_sum_array(tlen, treatment)
                ## Calculate psi score
                scores = []
                for i in xrange(1, tlen):
                    r_plus = treatment[i]
                    r_minus = control[i]
                    [ws, wr_plus, wr_minus] = determine_ws(i, csum, tsum)
                    scores.append(ws * (r_plus - r_minus) / (wr_plus + wr_minus) if wr_plus + wr_minus > 0.0 else 0.0)
                scores.append(0.0)
                # Output results
                fout.write(" ".join(["{:.6g}".format(x) for x in scores]) + "\n")

                if line_no % 1000 == 0:
                    print("FIN {}.".format(line_no))
