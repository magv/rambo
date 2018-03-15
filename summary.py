#!/usr/bin/env python3

import collections
import glob
import math
import sys

Integral = collections.namedtuple("Integral", "name d n alpha value error prob")

def nfmt(mean, sdev):
    sdev_digits = int(math.floor(math.log10(sdev)))
    sdev_digit0 = "%.0f" % (sdev*(1.0 + sys.float_info.epsilon)*10**(-sdev_digits))
    sdev_digits += len(sdev_digit0) - 1
    sdev_digit0 = sdev_digit0[0]
    mean_str = "%.0f" % (mean*(1.0 + sys.float_info.epsilon)*10**(-sdev_digits))
    mean_digits = len(mean_str) + sdev_digits
    if sdev_digits >= 0:
        return "%s(%s)%s" % (mean_str, sdev_digit0, "0"*sdev_digits)
    if mean_digits > 0:
        return "%s.%s(%s)" % (mean_str[:sdev_digits], mean_str[sdev_digits:], sdev_digit0)
    return "0.%s%s(%s)" % ("0"*(-mean_digits), mean_str, sdev_digit0)

def printf(fmt, *args):
    print(fmt.format(*args), end="")

integrals = {}

for filename in glob.glob("results/*.txt"):
    i_name = None
    i_d = None
    i_n = None
    i_value = None
    i_error = None
    i_alpha = None
    i_prob = None
    with open(filename, "r") as f:
        for line in f:
            if " = " in line:
                key, value = line.split(" = ", 1)
                key, value = key.strip(), value.strip()
                if key == "INTEGRAL": i_name = value
                if key == "D": i_d = value
                if key == "N": i_n = value
                if key == "integral": i_value = value
                if key == "error": i_error = value
                if key == "alpha": i_alpha = value
                if key == "prob": i_prob = value
    if i_alpha != "2^(-100)": continue
    if float(i_prob) > 0.9: continue
    integrals[i_name, i_d] = Integral(i_name, i_d, i_n, i_alpha, i_value, i_error, i_prob)

names = sorted(set(name for name, dim in integrals.keys()), key=lambda N: int(N))
dims = sorted(set(dim for name, dim in integrals.keys()))

printf("{:>3}", "")
for dim in dims:
    printf("{:>12}", "D=" + dim)
printf("\n")

for name in names:
    printf("{:>3}", name)
    for dim in dims:
        i = integrals.get((name, dim), None)
        if i is not None:
            printf("{:>12}", nfmt(float(i.value.replace("(", "").replace(")", "")), float(i.error)))
        else:
            printf("{:>12}", "--")
    printf("\n")
