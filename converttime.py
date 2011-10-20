#! /usr/bin/env python

times = []
t = raw_input()
while True:
    if len(t) > 0:
        times.append(t.strip())
    try:
        t = raw_input()
    except EOFError:
        break

times2 = [t.split('m') for t in times]
minutos = [int(t[0]) for t in times2]
segundos = [float(t[1][:-1]) for t in times2]
times = [t[0]*60 + t[1] for t in zip(minutos, segundos)]

for t in times:
    print "%.02f" % t

