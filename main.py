import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from PIL import Image
import numpy as np

def find_nearest(array, value):
    array = np.asarray(array)
    idx = 0
    minv = np.abs(array[0] - value)
    for i in range(1, len(array)):
        if np.abs(array[i] - value) < minv:
            idx = i
            minv = np.abs(array[i] - value)

    return idx

def getXYpos(relativelat, relativelon, platitude, plongitude):
    """ Calculates X and Y distances in meters.
    """
    deltaLatitude = platitude - relativelat
    deltaLongitude = plongitude - relativelon
    latitudeCircumference = 40075160 * np.cos(np.radians(relativelat))
    resultX = deltaLongitude * latitudeCircumference / 360
    resultY = deltaLatitude * 40008000 / 360
    return resultX, resultY

timestamps = []
xcoords = []
ycoords = []
R = 6371000

positionFile = open("position.txt", "r")
content = positionFile.readlines()
counter = 0
for line in content:
    if len(line) < 33:
        break
    timestamps.append(float(line[0:9]))
    ycoords.append(float(line[11:21]))
    xcoords.append(float(line[23:34]))

for i in range(len(timestamps)):
    temp = timestamps[i] - int(timestamps[i])
    temp += (int(timestamps[i]) % 100)
    temp += ((int(timestamps[i]) // 100) % 100) * 60
    temp += (int(timestamps[i]) // 10000) * 3600
    timestamps[i] = temp

startingPosition = [xcoords[0], ycoords[0], timestamps[0]]
for i in range(len(timestamps)):
    xcoords[i], ycoords[i] = getXYpos(startingPosition[1] / 60., startingPosition[0] / 60.,
                                      ycoords[i] / 60., xcoords[i] / 60.)
    xcoords[i] = 0 - xcoords[i]
    timestamps[i] -= startingPosition[2]

displacement = np.copy(timestamps)
for i in range(len(timestamps)):
    displacement[i] = np.sqrt((xcoords[0] - xcoords[i]) ** 2 +
                              (ycoords[0] - ycoords[i]) ** 2)

distance = np.copy(timestamps)
distance[0] = 0.
for i in range(1, len(timestamps)):
    distance[i] = distance[i-1] + np.sqrt((xcoords[i] - xcoords[i - 1]) ** 2 +
                                          (ycoords[i] - ycoords[i - 1]) ** 2)

# timestamps = [0., 0.1, 0.2, 0.3, 0.4]
# xcoords = [0, 1, 2, -1, -3]
# ycoords = [0, 1, -2, -6, -11]

root = tk.Tk()

verts = np.zeros((len(xcoords), 2))
for i in range(len(xcoords)):
    verts[i][0] = xcoords[i]
    verts[i][1] = ycoords[i]

fig = plt.figure()
ax = fig.add_subplot(111)
patch = ptch.Polygon(verts, closed=False, lw=2, fc='none', ec="black")
ax.add_patch(patch)
ax.set_xlim(np.min(xcoords) - (np.max(xcoords) - np.min(xcoords))/50, np.max(xcoords) + (np.max(xcoords) - np.min(xcoords))/50)
ax.set_ylim(np.min(ycoords) - (np.max(ycoords) - np.min(ycoords))/50, np.max(ycoords) + (np.max(ycoords) - np.min(ycoords))/50)
plt.plot([xcoords[0]], [ycoords[0]], marker='o', markersize=5)
fig.savefig("temp1.png")
im = Image.open("temp1.png")
im.save("temp1.ppm")
mainGraph = tk.PhotoImage(file="temp1.ppm")
plt.close()
time = 50.


# set weights for each row and column
root.columnconfigure(0, pad=3, weight=1)
root.rowconfigure(0, pad=3, weight=0)
root.rowconfigure(1, pad=3, weight=40)
root.rowconfigure(2, pad=3, weight=0)

# row 1 is taken up by the title
title = tk.Label(root, text="Probe Application")
title.grid(row=0, sticky=tk.N+tk.E+tk.S+tk.W)

# graphFrame contains 2 images (temporarily buttons)
graphFrame = tk.Frame(root)
graphFrame.rowconfigure(0, pad=0, weight=1)
graphFrame.columnconfigure(0, pad=2, weight=1)
graphFrame.columnconfigure(1, pad=2, weight=1)
graphFrame.grid(row=1, sticky=tk.N+tk.E+tk.S+tk.W)

# image 1 is the path, a graph showing where the probe went

# image1 = tk.Canvas(graphFrame, width=100, height=100)
# image1.create_image(0, 0, image=mainGraph, anchor=tk.NW)
image1 = tk.Label(graphFrame, image=mainGraph)
image1.grid(row=0, column=0, sticky=tk.N+tk.E+tk.S+tk.W)

# image 2 is the additional info, showing measurements taken at each point in the path
tempButton2 = tk.Button(graphFrame, text="measurement graph", bg="green", height=30)
tempButton2.grid(row=0, column=1, sticky=tk.N+tk.E+tk.S+tk.W)

# sliderFrame lets the user decide where along the time axis to go
sliderValue = tk.DoubleVar()
slider = tk.Scale(root, from_=0., to=timestamps[len(timestamps) - 1], label="Time", length=1500, orient=tk.HORIZONTAL,
                  relief=tk.SOLID, resolution=0.01, sliderlength=20, variable=sliderValue, width=30)
slider.grid(row=2, sticky=tk.N+tk.E+tk.S+tk.W)

# textFrame holds two text boxes. One for general information and another for time/location specific information
textFrame = tk.Frame(root)
textFrame.rowconfigure(0, pad=0, weight=1)
textFrame.columnconfigure(0, pad=2, weight=1)
textFrame.columnconfigure(1, pad=2, weight=1)
textFrame.grid(row=3, sticky=tk.N+tk.E+tk.S+tk.W)

# text box 1 holds the time-specific info
specBox = tk.Text(textFrame)
specBox.insert(tk.END, "Time-specific info changed by the slider")
specBox.grid(row=0, sticky=tk.N+tk.E+tk.S+tk.W)

# text box 2 holds general info
genBox = tk.Text(textFrame)
genBox.insert(tk.END, "Distance Travelled : " + str(np.round(distance[len(timestamps) - 1], 2)) + " m")
genBox.insert(tk.END, "\nTotal Displacement : " + str(np.round(displacement[len(timestamps) - 1], 2)) + " m")
genBox.insert(tk.END, "\nTime Elapsed       : " + str(np.round(timestamps[len(timestamps) - 1], 2)) + " s")
genBox.grid(row=0, column=1, sticky=tk.N+tk.E+tk.S+tk.W)


def updateScale(sliderv):
    specBox.delete('1.0', tk.END)

    _fig = plt.figure()
    _ax = _fig.add_subplot(111)
    _patch = ptch.Polygon(verts, closed=False, lw=2, fc='none', ec="black")
    _ax.add_patch(_patch)
    _ax.set_xlim(np.min(xcoords) - (np.max(xcoords) - np.min(xcoords)) / 50,
                 np.max(xcoords) + (np.max(xcoords) - np.min(xcoords)) / 50)
    _ax.set_ylim(np.min(ycoords) - (np.max(ycoords) - np.min(ycoords)) / 50,
                 np.max(ycoords) + (np.max(ycoords) - np.min(ycoords)) / 50)

    # find the current selected point
    index = find_nearest(timestamps, float(sliderv))
    plt.plot([xcoords[index]], [ycoords[index]], marker='o', markersize=5)

    _fig.savefig("temp1.png")
    _im = Image.open("temp1.png")
    _im.save("temp1.ppm")
    global mainGraph
    mainGraph = tk.PhotoImage(file="temp1.ppm")
    image1.configure(image=mainGraph)
    plt.close()

    specBox.insert(tk.END, "Distance Travelled so far : " + str(np.round(distance[index], 2)) + " m")
    specBox.insert(tk.END, "\nDisplacement from start   : " + str(np.round(displacement[index], 2)) + " m")


slider.configure(command=updateScale)


root.mainloop()


