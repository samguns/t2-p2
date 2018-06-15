import csv
import numpy as np
import matplotlib.pyplot as plt

radar_nis_array = []
laser_nis_array = []

with open('./nis.csv') as csvfile:
    reader = csv.reader(csvfile)
    for line in reader:
        if line[0] is not '':
            radar_nis_array.append(line[0])

        if line[1] is not '':
            laser_nis_array.append(line[1])


radar_nis = np.array(radar_nis_array, dtype=float)
laser_nis = np.array(radar_nis_array, dtype=float)

radar_above_95_thresh = np.where(radar_nis > 7.82)
laser_above_95_thresh = np.where(laser_nis > 5.99)

radar_above_rate = radar_above_95_thresh[0].size / radar_nis.size
laser_above_rate = laser_above_95_thresh[0].size / laser_nis.size

plt.plot(radar_nis)
plt.axhline(7.82, color='red')
plt.title("Radar NIS - " + str(radar_above_rate*100) + "% above 7.82")
plt.show()

plt.plot(laser_nis)
plt.axhline(5.99, color='red')
plt.title("Laser NIS - " + str(laser_above_rate*100) + "% above 5.99")
plt.show()