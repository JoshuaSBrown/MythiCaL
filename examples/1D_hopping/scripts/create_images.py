import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat

site=[]
dwell_Time=[]

with open("Trajectory.txt") as fid:

    line = fid.readline()
    words = line.split()
    num_sites = int(words[3])

    line = fid.readline()
    words = line.split()
    iterations = int(words[3])

    line = fid.readline()
    line = fid.readline()
    

    for line in fid:
        site.append(int(line.split()[0]))
        dwell_Time.append(float(line.split()[1]))

sites = np.asarray(site)

for i in range(0,iterations):
    fig = plt.figure(figsize = (num_sites+1,5))
    ax = plt.gca()
    for j in range(1,num_sites+1):
        draw_circ = plt.Circle((j,0),0.4, color='r')
        ax.add_artist(draw_circ)
      
    draw_circ = plt.Circle((sites[i],0),0.4, color='b')
    ax.add_artist(draw_circ)

    ax.set_xlim(0, num_sites+1)
    ax.set_ylim(-2,2)
    ax.set_aspect('equal')
    plt.xlabel("x axis")
    plt.savefig("hop_" + str(i) + ".png")

    plt.close()
