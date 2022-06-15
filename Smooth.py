# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 11:50:39 2018

@author: Walid.Hanifi
"""

# Libraries
from datetime import datetime
startTime = datetime.now()
import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.interpolate as sp
import matplotlib.patches as mpatches

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
INPUT
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# COMPULSORY SETTINGS
STRAIGHT_LINE_SECTIONS = 10
BENDING_RADIUS = 30

# OPTIONAL SETTINGS
major_bends1 = 0.001             #default = 0.1
major_bends2 = 0                #default = 0
plot_circles = True             #default = True

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
SMOOTHING ALGORITHM
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Load variables used in algorithm.py from variable pass.txt
f = open("variable pass.txt", "r")

Output_File_Name_Path = f.readline().rstrip('\n')
cell_size = float(f.readline().rstrip('\n'))
print(cell_size)
min_x = float(f.readline().rstrip('\n'))
min_y = float(f.readline().rstrip('\n'))
max_x = float(f.readline().rstrip('\n'))
max_y = float(f.readline().rstrip('\n'))
delta_start_x = float(f.readline().rstrip('\n'))
delta_start_y = float(f.readline().rstrip('\n'))
delta_end_x = float(f.readline().rstrip('\n'))
delta_end_y = float(f.readline().rstrip('\n'))
f.close()

## MISCELLANEOUS CHECKS
# Closes any open figures, to replot the smooth figure
plt.close('all')
fig, ax = plt.subplots(figsize=(12, 8))
try:
    ax.imshow(img3, extent=(0, NoOfColumns, 0, NoOfRows), interpolation = 'none', filternorm = 1, cmap=colormap)
except NameError:
    sys.exit("RUN Algorithm.py FIRST!")
plt.xticks([])
plt.yticks([])
plt.title("Path" + str(Output_File_Name_Path), fontsize=7)
# Calibrating the given Bending radius to match our graphs coordinates
RADIUS = BENDING_RADIUS / cell_size
# Checks whether any clicked points were saved or not, and whether to use it or not, if saved_coords text file is empty then smooth according to the specified number of splines
if (os.stat('saved_coords.txt').st_size==0) == True:
    saved_coords = np.loadtxt("Path" + str(Output_File_Name_Path) + ".xyz", delimiter=',')
    saved_coords[:,0] = (saved_coords[:,0] - min_x)/cell_size
    saved_coords[:,1] = (saved_coords[:,1] - min_y)/cell_size
# If saved_coords text file is not empty then smooth according to the coordinates in the saved_coords text file
elif (os.stat('saved_coords.txt').st_size==0) == False:
    saved_coords = np.loadtxt('saved_coords.txt')
x = saved_coords[:,0]
y = saved_coords[:,1]

## PERFORM THE SMOOTHING
# Using splines to smooth the jagged path, into a set of straight line sections
new_y = np.linspace(y.min(), y.max(), (STRAIGHT_LINE_SECTIONS + 1))
new_x = sp.interp1d(y, x, kind='linear')(new_y)
old_x = new_x   #initialising for the bending algorithm
old_y = new_y   #initialising for the bending algorithm



"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
BENDING FUNCTION INITIALISATION
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def bending(point1x, point1y, point2x, point2y, point3x, point3y):

    ## INITIALISING VARIABLES
    LINE1_X = point1x
    LINE1_Y = point1y
    INTERSECTION_X = point2x
    INTERSECTION_Y = point2y
    LINE2_X = point3x
    LINE2_Y = point3y
    POINT1 = LINE1_X, LINE1_Y
    POINT2 = INTERSECTION_X, INTERSECTION_Y
    POINT3 = LINE2_X, LINE2_Y
    LINE1 = ((LINE1_X, LINE1_Y),(INTERSECTION_X, INTERSECTION_Y))
    LINE2 = ((LINE2_X, LINE2_Y),(INTERSECTION_X, INTERSECTION_Y))
    
    ## AUTOMATED CIRCLE PLACEMENT
    # Automatic detection of which side to plot the circle on, using the atan2 function
    Dir_C_to_A = math.atan2(point1y - point2y, point1x - point2x)
    Dir_C_to_B = math.atan2(point3y - point2y, point3x - point2x)
    Angle_ACB = Dir_C_to_A - Dir_C_to_B
    Pi = math.pi
    if (Angle_ACB > Pi):
        Angle_ACB -= 2*Pi
    elif (Angle_ACB < -Pi):
        Angle_ACB += 2*Pi
    # Converting the Angle to degrees, and determining the sign
    if (np.degrees(Angle_ACB)) > 0:
        reverse = True
    elif (np.degrees(Angle_ACB)) < 0:
        reverse = False
    # Initialising Functions for the Angle between two lines and two points
    def dot(vA, vB):
        return vA[0]*vB[0]+vA[1]*vB[1]
    def ang_lines(lineA, lineB):
        # Get vector form
        vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
        vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
        # Get dot product
        dot_prod = dot(vA, vB)
        # Get magnitudes
        magA = dot(vA, vA)**0.5
        magB = dot(vB, vB)**0.5
        # Get angle in radians and then convert to degrees
        try:
            angle = math.acos(dot_prod/magB/magA)
        except ValueError:
            angle = 0
        ang_deg = math.degrees(angle)%360
        if ang_deg - 180 >= 0:
            return 360 - ang_deg
        else: 
            return ang_deg
    def ang_points(pointA, pointB):
        dx = abs(pointA[0]-pointB[0])
        dy = abs(pointA[1]-pointB[1])
        ang_new = math.degrees(math.atan2(dy,dx))
        return ang_new
    # Miscellanous angle calculations
    beta = (ang_lines(LINE1, LINE2))/2
    theta = 180 - ang_lines(LINE1, LINE2)
    gamma = (180 - (beta * 2))/2
    l1_length = RADIUS/(np.cos(np.radians(gamma)))   # length from centre of circle to point 2
    YAXIS = ((0,0),(0,1000)) # Initialising a yaxis line
    default_theta = ((180 - ang_lines(YAXIS, LINE2)) - ang_lines(YAXIS, LINE1))
    # Shape 1:
    if default_theta > 0:
        yaxis_line1 = - ang_lines(YAXIS, LINE1)
        E = ang_lines(YAXIS, LINE1)
        sign = 1
    # Shape 2:  
    elif default_theta < 0 and (abs(round(default_theta,1)) == abs(round(theta,1))):
        yaxis_line1 = 0
        E = 0
        sign = -1
    # Shape 3:
    elif default_theta < 0 and (abs(round(default_theta,1)) != abs(round(theta,1))):
        yaxis_line1 = 0
        E = 0
        sign = 1
    theta = (180 - ang_lines(YAXIS, LINE2)) + yaxis_line1
    l1_angle = E + sign * theta + beta
    # Choosing which side  of the path to plot the circle on
    if reverse == False:
        l1_angle_rads = np.radians(-l1_angle)
    elif reverse == True:
        l1_angle_rads = np.radians(l1_angle)
    POINT4 = INTERSECTION_X, INTERSECTION_Y + l1_length
    def rotate(origin, point, angle):
        #Rotate a point counterclockwise by a given angle around a given origin. The angle should be given in radians.
        ox, oy = origin
        px, py = point
        qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
        qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
        return qx, qy
    CIRCLE_CENTER = rotate(POINT2, POINT4, l1_angle_rads)
    # Calculating and then Plotting the arcs
    new_extended_x = CIRCLE_CENTER[0] + RADIUS
    x_ref_line = ((CIRCLE_CENTER[0], CIRCLE_CENTER[1]),(new_extended_x, CIRCLE_CENTER[1]))
    first_line = ((CIRCLE_CENTER[0], CIRCLE_CENTER[1]),(INTERSECTION_X, INTERSECTION_Y))
    the_angle = ang_lines(x_ref_line, first_line)
    # Default Library for all circles
    if reverse == False:
        new_angle = -gamma + the_angle
    elif reverse == True:
        new_angle = -gamma - the_angle
    # Developed Library for specific cases
    if reverse == True and point3x > point2x > point1x:
        new_angle = -gamma - the_angle
    elif reverse == True and point3x < point2x and point3x < point1x:
        new_angle = -gamma + the_angle
    if reverse == False and point3x < point2x < point1x:
        new_angle = -gamma - the_angle

    ## ITERATION 2 INITIALISATION
    # Making the script compatible with a second iteration
    Arc_Length = abs(2*np.radians(gamma)*RADIUS)
    Circumference = 2*RADIUS*np.pi    
    if Arc_Length > major_bends*Circumference:
        pac = mpatches.Arc((CIRCLE_CENTER), RADIUS*2, RADIUS*2, angle=(new_angle), theta1=0, theta2=2*gamma, edgecolor=color_circle)
        plt.gcf().gca().add_artist(pac)
        # Append the corner point to an array, so we can join them up later
        corners_x.append(point2x)
        corners_y.append(point2y)
        # Plot the full circles or just the arc?
        if plot_circles == True:
            circle1 = plt.Circle((CIRCLE_CENTER), radius=RADIUS, fill=False, color=color_circle)
            plt.gcf().gca().add_artist(circle1)
        elif plot_circles == False:
            0
    if Arc_Length < major_bends*Circumference:
        if plot_circles == True:
            circle1 = plt.Circle((CIRCLE_CENTER), radius=RADIUS, fill=False, color=color_circle)
            plt.gcf().gca().add_artist(circle1)
        elif plot_circles == False:
            0
    return corners_x, corners_y

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
LOOPING THE BENDS
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## 1ST ITERATION
# Running the Bending Algorithm for the 1st iteration
major_bends = major_bends1
color_path = '#0080ff'
color_circle = '#99c2ff'
corners_x = []  #Initialising corner arrays (for iteration2)
corners_y = []
# Attaching the smoothed path to the start/end point correctly
old_x[0] = delta_start_x
old_y[0] = delta_start_y
old_x = old_x[::-1]     #temporarily flipping old_x
old_x[0] = delta_end_x
old_x = old_x[::-1]
old_y = old_y[::-1]     #temporarily flipping old_y
old_y[0] = delta_end_y
old_y = old_y[::-1]
# Plotting the path
plt.plot(new_x, new_y, color=color_path)
for i in range(len(old_x)):
    if (i+1) >= len(old_x) or (i+2) >= len(old_x):
        break
    elif (i+1) < len(old_x) or (i+2) < len(old_x):
        point1x = old_x[i]
        point2x = old_x[i+1]
        point3x = old_x[i+2]
        point1y = old_y[i]
        point2y = old_y[i+1]
        point3y = old_y[i+2]
        bending(point1x, point1y, point2x, point2y, point3x, point3y)

## 2ND ITERATION
# Running the Bending Algorithm for the 2nd iteration
major_bends = major_bends2
color_path = '#ff0000'
color_circle = '#ff0000'
corners = np.column_stack((corners_x, corners_y))
old_x2 = corners[:,0]
old_y2 = corners[:,1]
# Fixing Start point being cutoff
first_x = old_x[0]
first_y = old_y[0]
old_x2 = np.insert(old_x2, 0, first_x)
old_y2 = np.insert(old_y2, 0, first_y)
# End point
last_x = old_x[::-1]
last_x = last_x[0]
last_y = old_y[::-1]
last_y = last_y[0]
old_x2 = np.append(old_x2, last_x)
old_y2 = np.append(old_y2, last_y)
# Plotting the path
plt.plot(old_x2, old_y2, color=color_path)
for i in range(len(old_x2)):
    if (i+1) >= len(old_x2) or (i+2) >= len(old_x2):
        break
    elif (i+1) < len(old_x2) or (i+2) < len(old_x2):
        point1x = old_x2[i]
        point2x = old_x2[i+1]
        point3x = old_x2[i+2]
        point1y = old_y2[i]
        point2y = old_y2[i+1]
        point3y = old_y2[i+2]
        bending(point1x, point1y, point2x, point2y, point3x, point3y)

## MISCELLANEOUS ITEMS (LEGEND, START/END CIRCLES)
# Plotting the start and end points
circle_start = plt.Circle((delta_start_x, delta_start_y), radius=1, fill=True, color='#00ff00')
plt.gcf().gca().add_artist(circle_start)
circle_end = plt.Circle((delta_end_x, delta_end_y), radius=1, fill=True, color='#ffff00')
plt.gcf().gca().add_artist(circle_end)
# Placing the legend
blue_iter1 = mpatches.Patch(color='#0080ff', label='Iteration 1')
red_iter2 = mpatches.Patch(color='#ff0000', label='Iteration 2')
plt.legend(handles=[blue_iter1, red_iter2])

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
SAVING SMOOTHED PATH
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## COORDINATE CALIBRATION
new_x = new_x * cell_size
new_y = new_y * cell_size
new_x = new_x + min_x
new_y = new_y + min_y

## OUTPUTTING .XYZ PATH
# Creating an empty Z column and attaching it to our xy coordinates
length = len(new_x)
ones = list(np.random.randint(1, 3, size=length))
ones = np.asarray(ones)
ones = np.reshape(ones, (length, 1))
new_xy = np.column_stack((new_x, new_y))
new_xyz = np.column_stack((new_xy, ones))
# Saving the Smoothed path in an XYZ coordinates file
np.savetxt("Path_Smoothed" + str(Output_File_Name_Path) + ".xyz", new_xyz, fmt='%f', delimiter=',')

## MISCELLANEOUS
# Set the output plot to be equal in size
plt.gca().set_aspect('equal', 'datalim')

print((datetime.now() - startTime), ': TOTAL TIME')