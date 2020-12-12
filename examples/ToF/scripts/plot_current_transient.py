
import numpy as np
import matplotlib.pyplot as plt

moving_avg = 5

with open("transient_current.txt",'r') as tc:
    tc.readline() # Ignore header
    num_pts = int(tc.readline())
    time = np.zeros(num_pts)
    current = np.zeros(num_pts)
    current_moving_avg = np.zeros(num_pts)
    charges = np.zeros(num_pts)

    lines = tc.readlines()
    count = 0
    # Strips the newline character 
    for line in lines: 
        words = line.split()
        time[count] = float(words[0])
        current[count] = float(words[1])
        charges[count] = int(words[2])
        count = count + 1

    value = 0.0
    for ind in range(0,num_pts):
        value = value + current[ind]
        if ind < moving_avg:
            current_moving_avg[ind] = value/float(ind+1)
        elif ind > num_pts - moving_avg:
            value = value - current[ind-moving_avg]
            current_moving_avg[ind] = value/float(num_pts - ind)
        else:
            value = value - current[ind-moving_avg]
            current_moving_avg[ind] = value/float(moving_avg)

    xlim_val = 2E-9
    plt.plot(time,current)
    plt.plot(time,current_moving_avg)
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    plt.xlabel("Time $[s]$")
    plt.ylabel("Displacement Current $[A / cm^2]$")
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig("DisplacementCurrentLogLog.png")
    plt.xscale('linear')
    plt.yscale('linear')
    ax = plt.gca()
    plt.xlim(0,xlim_val)
#    plt.ylim(0,20000)
    ax2 = ax.twinx()
    ax2.plot(time,charges,color='g')
    ax2.set_ylabel("Number of Charges Remaining")
    plt.ylim(0,210)
    plt.savefig("DisplacementCurrent.png")
