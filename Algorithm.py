# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:03:55 2018

@author: Walid.Hanifi
"""

# Libraries
print('CALCULATION TIME:')
from datetime import datetime
startTime = datetime.now()
import sys
import os
import os.path
import math
import numpy as np
from numpy import median
import matplotlib.pyplot as plt
import pandas as pd   # Must be Installed
import networkx as nx   # Must be Installed
from tqdm import tqdm   # Must be Installed

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
INPUT
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# COMPULSORY SETTINGS
Working_Folder = r'C:\Users\WHanifi\Dropbox\2. Programming\1. Automated Routing'
Working_Folder = os.getcwd()
Project_Title = "test Project"
Filename = "test_input_raster.xyz"
START_X, START_Y = [671572.00, 5173390.00]
END_X, END_Y = [671892.00, 5173738.00]
Delimiter = " "
Heuristic = "Dijkstra"
        # Heuristic Options:
                            #A_Star
                            #Dijkstra
                            #Bidirectional_Search
                            #Bellman_Ford

colormap = 'ocean'        #See library website for all options: https://matplotlib.org/stable/tutorials/colors/colormaps.html

# OPTIONAL SETTINGS
Directional_Weighting = 0.0005          # controls the behaviour behind selecting the nearest edge, default = 0.001
Elevation_Weighting = 1                 # default = 1
number_of_derivatives = 1               # 0 = map based on z elevation, 1 = map based on slopes
dpi = 1000                              # image output dpi, default = 1000
Special_Z_elevation_Map = True          # default = True (uses a fancy z elevation map)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
CONVERTING XYZ FILE TO WORKABLE RASTER
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Change the Working Folder to the one specified
os.chdir(Working_Folder)

# Delete the 'saved_coords.txt' file if it exists
open('saved_coords.txt', 'w').close()

# Close all open plots
plt.close('all')

# Seperating xyz file into x, y and z columns
file = pd.read_csv(Filename, delimiter=Delimiter, header=None)
file = file.values
x_coords = file[:,0]
y_coords = file[:,1]
z_coords = file[:,2]


# Finding the bottom left and top right coordinates
min_x = min(x_coords) #left
max_x = max(x_coords) #right
min_y = min(y_coords) #bottom
max_y = max(y_coords) #top

# Finding Cell size from x column differences
cell_size = abs(x_coords[1] - x_coords[0])

# Readjusting coordinates to start from 0,0
new_max_x = (max_x - min_x)/cell_size
new_max_y = (max_y - min_y)/cell_size
new_max_x = new_max_x.astype(int)
new_max_y = new_max_y.astype(int)

# Defining number of rows and columns of grid from this xyz
rows = new_max_y
columns = new_max_x

# Creating a new xyz, starting from 0,0, with intervals of 1
x_coords_new = (file[:,0]-min_x)/cell_size
y_coords_new = (file[:,1]-min_y)/cell_size
new_file = np.column_stack((x_coords_new, y_coords_new, z_coords))

# Converting the new xyz into a raster
df = pd.DataFrame(new_file, columns=['x', 'y', 'z'])
df = df.pivot_table(values='z', index='x', columns='y')
num = df.values

# Rotating this array to be properly shown
num = np.rot90(num)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
EDGE COST CALCULATION
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
startTime0 = datetime.now()

    # INITIALISING OPTIONS
fig, ax = plt.subplots(figsize=(12, 8))
NoOfRows, NoOfColumns = num.shape   # Grabbing the NoOfRows and NoOfColumns from our bathymetry

# Define the Missing Data and ymax variables
num = np.nan_to_num(num)            #Converts all missing data (nan) to zeros
Set_Missing_Data = median(num)      #uses the median value of the raster as the missing data value so that colormap stays consistent
ymax = Set_Missing_Data + 30        #the maximum y value of the on bottom roughness plot

# Set the missing data as a variable, from the weighted raster
WeightedRasterWithoutBorder = num
WeightedRasterWithoutBorder[WeightedRasterWithoutBorder == 0] = Set_Missing_Data

# Reshaping the loaded array to work with our calculations
WeightedRasterWithoutBorder = np.reshape(WeightedRasterWithoutBorder, (NoOfRows, NoOfColumns))     

# Creating a Nodal Grid (Top left, to bottom right corner) to be used in later calculations
NumberOfNodes = (NoOfRows * NoOfColumns) + 1
NodalRasterWithoutBorder = range(1, NumberOfNodes)
NodalRasterWithoutBorder = list(NodalRasterWithoutBorder)
NodalRasterWithoutBorder = np.reshape(NodalRasterWithoutBorder, (NoOfRows, NoOfColumns))

# Setting aside an untampered weighted raster
WeightedRasterUntampered = WeightedRasterWithoutBorder

# Calculating the slope of the weighted raster
if number_of_derivatives > 0:
    WeightedRasterWithoutBorder = np.diff(WeightedRasterWithoutBorder, number_of_derivatives)
    WeightedRasterWithoutBorder = np.hstack((WeightedRasterWithoutBorder, np.tile(WeightedRasterWithoutBorder[:, [-1]], 2)))
    WeightedRasterWithoutBorder = WeightedRasterWithoutBorder[:,:(len(WeightedRasterWithoutBorder[0]) -(- number_of_derivatives + 2))] # Derivatives > 2 crash here
    slope_string = "slope"
elif number_of_derivatives == 0:
    slope_string = "z-elevation"
elif number_of_derivatives < 0:
    sys.exit("INCORRECT NUMBER OF DERIVATIVES SPECIFIED!")

# Adding borders to the Nodal and Weighted Rasters, so that the graph does not search over these borders
NodalRaster = np.pad(NodalRasterWithoutBorder, 1, 'constant', constant_values=0)
WeightedRaster = np.pad(WeightedRasterWithoutBorder, 1, 'constant', constant_values=999999)

# The number of rows and columns including the number of borders that were previously added
new_rows = WeightedRaster.shape[0]
new_cols = WeightedRaster.shape[1]

# Creating a factored version of the cell size, to later use when relaxing edge cost directional calculations
cell_size_factored = cell_size * Directional_Weighting

# Converting Start and End X, Y points to Node format
delta_start_x = (START_X - min_x)/cell_size
delta_start_y = (- START_Y + max_y)/cell_size
START_NODE = NodalRasterWithoutBorder[delta_start_y.astype(int), delta_start_x.astype(int)] #TOO LARGE/TOO SMALL GRIDS CRASH HERE
delta_end_x = (END_X - min_x)/cell_size
delta_end_y = (- END_Y + max_y)/cell_size
END_NODE = NodalRasterWithoutBorder[delta_end_y.astype(int), delta_end_x.astype(int)]

# Initialising the output file name
Output_File_Name_Path = " ["+str(Project_Title)+", "+str(Heuristic)+"] (dydx="+str(number_of_derivatives)+", dw="+str(Directional_Weighting)+", ew="+str(Elevation_Weighting)+", cells="+str(round(cell_size))+"m)"
Output_File_Name_Edge = " ["+str(Project_Title)+"] (dydx="+str(number_of_derivatives)+", dw="+str(Directional_Weighting)+", ew="+str(Elevation_Weighting)+", cells="+str(round(cell_size))+"m)"

# Calculating the Edge Costs for each node, in all eight directions and stores these into EdgeCosts.txt
if os.path.isfile("EdgeCosts" + str(Output_File_Name_Edge) + ".txt") == False:
    with open("EdgeCosts" + str(Output_File_Name_Edge) + ".txt", "w") as text_file:
        for i in tqdm(range(new_rows), ncols=100, desc="Edge Costs", unit='B', unit_scale=True, unit_divisor=1024):
            for j in range(new_cols):
                initial_node_number = NodalRaster[i,j]
                initial_weighting = WeightedRaster[i,j]
                try:
                    offset_node_number_one = NodalRaster[i-1,j-1]
                    offset_node_number_two = NodalRaster[i-1,j]
                    offset_node_number_three = NodalRaster[i-1,j+1]
                    offset_node_number_four = NodalRaster[i,j+1]
                    offset_node_number_five = NodalRaster[i+1,j+1]
                    offset_node_number_six = NodalRaster[i+1,j]
                    offset_node_number_seven = NodalRaster[i+1,j-1]
                    offset_node_number_eight = NodalRaster[i,j-1]
                except IndexError:
                    0
                try:
                    offset_weighting_one = WeightedRaster[i-1,j-1]
                    offset_weighting_two = WeightedRaster[i-1,j]
                    offset_weighting_three = WeightedRaster[i-1,j+1]
                    offset_weighting_four = WeightedRaster[i,j+1]
                    offset_weighting_five = WeightedRaster[i+1,j+1]
                    offset_weighting_six = WeightedRaster[i+1,j]
                    offset_weighting_seven = WeightedRaster[i+1,j-1]
                    offset_weighting_eight = WeightedRaster[i,j-1]
                except IndexError:
                    0
    #The Formula for calculating the Edge Weightings in all eight directions
                edge_weighting_one = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_one))**2 + cell_size_factored**2 + cell_size_factored**2)))
                edge_weighting_two = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_two))**2 + cell_size_factored**2)))
                edge_weighting_three = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_three))**2 + cell_size_factored**2 + cell_size_factored**2)))
                edge_weighting_four = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_four))**2 + cell_size_factored**2)))
                edge_weighting_five = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_five))**2 + cell_size_factored**2 + cell_size_factored**2)))
                edge_weighting_six = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_six))**2 + cell_size_factored**2)))
                edge_weighting_seven = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_seven))**2 + cell_size_factored**2 + cell_size_factored**2)))
                edge_weighting_eight = math.sqrt(abs(((Elevation_Weighting*(initial_weighting - offset_weighting_eight))**2 + cell_size_factored**2)))
                
    #Saving these edge cost values to a string in the EdgeCosts.txt file
                edge_one_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_one) + ", " +  str(edge_weighting_one) + "),"
                edge_two_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_two) + ", " +  str(edge_weighting_two) + "),"
                edge_three_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_three) + ", " +  str(edge_weighting_three) + "),"
                edge_four_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_four) + ", " +  str(edge_weighting_four) + "),"
                edge_five_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_five) + ", " +  str(edge_weighting_five) + "),"
                edge_six_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_six) + ", " +  str(edge_weighting_six) + "),"
                edge_seven_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_seven) + ", " +  str(edge_weighting_seven) + "),"
                edge_eight_string = "(" + str(initial_node_number) + ", " + str(offset_node_number_eight) + ", " +  str(edge_weighting_eight) + "),"
               
                if edge_weighting_one < 900 and initial_node_number != 0 and offset_node_number_one != 0:
                    print(edge_one_string, file=text_file)
                if edge_weighting_two < 900 and initial_node_number != 0 and offset_node_number_two != 0:
                    print(edge_two_string, file=text_file)
                if edge_weighting_three < 900 and initial_node_number != 0 and offset_node_number_three != 0:
                    print(edge_three_string, file=text_file)
                if edge_weighting_four < 900 and initial_node_number != 0 and offset_node_number_four != 0:
                    print(edge_four_string, file=text_file)
                if edge_weighting_five < 900 and initial_node_number != 0 and offset_node_number_five != 0:
                    print(edge_five_string, file=text_file)
                if edge_weighting_six < 900 and initial_node_number != 0 and offset_node_number_six != 0:
                    print(edge_six_string, file=text_file)
                if edge_weighting_seven < 900 and initial_node_number != 0 and offset_node_number_seven != 0:
                    print(edge_seven_string, file=text_file)
                if edge_weighting_eight < 900 and initial_node_number != 0 and offset_node_number_eight != 0:
                    print(edge_eight_string, file=text_file)
    
    # Open up our new EdgeCosts.txt file
    with open("EdgeCosts" + str(Output_File_Name_Edge) + ".txt", 'r') as text_file :
      filedata = text_file.read()
    
    # Correct the formatting
    filedata = filedata.replace('.0,', ',')
    filedata = filedata.replace('(', '')
    filedata = filedata.replace(', ', ' ')
    filedata = filedata.replace('),', '')
    
    # Save the file out again
    with open("EdgeCosts" + str(Output_File_Name_Edge) + ".txt", 'w') as text_file:
      text_file.write(filedata)
    networkx = pd.read_csv("EdgeCosts" + str(Output_File_Name_Edge) + ".txt", delimiter=' ', header=None)   #sep = delimiter for pandas
    networkx = networkx.values #converts the dataframe to a normal array

print((datetime.now() - startTime0), ': Edge Costs Time Taken')

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
PATHFINDING ALGORITHM (USING THE NETWORKX LIBRARY)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
startTime2 = datetime.now()
    
# Read Edgelist
G = nx.read_weighted_edgelist("EdgeCosts" + str(Output_File_Name_Edge) + ".txt", nodetype=float)

# Choose Heuristic, and run algorithm
if Heuristic=="Bellman_Ford":
    p = nx.bellman_ford_path(G, source=START_NODE, target=END_NODE)
    length = nx.bellman_ford_path_length(G, source=START_NODE, target=END_NODE, weight='weight')
elif Heuristic=="Dijkstra":
    p = nx.dijkstra_path(G, source=START_NODE, target=END_NODE, weight='weight')
    length = nx.dijkstra_path_length(G, source=START_NODE, target=END_NODE, weight='weight')
elif Heuristic=="Bidirectional_Search":
    length,path = nx.bidirectional_dijkstra(G, source=START_NODE, target=END_NODE, weight='weight')
    p = path
elif Heuristic=="A_Star":
    0
    p = nx.astar_path(G, source=START_NODE, target=END_NODE, weight='weight')
    length = nx.astar_path_length(G, source=START_NODE, target=END_NODE, weight='weight')
else:
    sys.exit("HEURISTIC NOT RECOGNIZED!")

print((datetime.now() - startTime2), ': ' + str(Heuristic) + ' Time Taken')

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
PLOTTING
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
startTime3 = datetime.now()

# Creating Image 1: The Bathymetry Map
if Special_Z_elevation_Map == True:
    img1 = WeightedRasterUntampered
    img1 = img1.astype(float)
elif Special_Z_elevation_Map == False:
    img1 = WeightedRasterWithoutBorder
    img1 = img1.astype(float)

# Creating Image 2: The OutputPath with blank cells represented as zeros
img2_part1 = NodalRasterWithoutBorder
img2_part2 = np.asarray(p)
img2 = np.in1d(img2_part1,img2_part2)    #Combining img2_part1 and img2_part2 by matching together
img2 = img2.astype(int)
img2 = np.reshape(img2, (NoOfRows, NoOfColumns))            
img2[img2 == 1] = 1000

# Combines Image 1 and Image 2 to get Image 3: The Output Path on the Bathymetry Map
img3 = np.where(img2==0, img1, img2) #Ignore the zeros when combining the images
img3 = np.ma.masked_where(img3==1000, img3)

# Plotting Image 3
ax.imshow(img3, extent=(0, NoOfColumns, 0, NoOfRows), interpolation = 'none', filternorm = 1, cmap=colormap)

# Return coordinates of point on map once clicked
ix = []
iy = []
coords = []
line, = ax.plot(ix, iy, marker="o")
def doubleclick(event):
    # Draws the button, saves its coordinates
    if event.button == 1 and event.dblclick:
        # Saves the clicked coordinates
        global ix, iy
        ix, iy = event.xdata, event.ydata
        print("xdata: " + str(event.xdata) + " ydata: " + str(event.ydata))
        global coords
        coords.append((ix, iy))
        np.savetxt('saved_coords.txt', coords, fmt='%f')
        # Creates a new button if the left click button is pressed
        line.set_xdata(ix)
        line.set_ydata(iy)
        fig.canvas.draw()
cid = fig.canvas.mpl_connect('button_press_event', doubleclick)

# Hide the x and y axis, to display a larger bathymetry map
plt.xticks([])
plt.yticks([])

# Display a title for the graph
plt.title("Path" + str(Output_File_Name_Path), fontsize=7)
plt.savefig("Path" + str(Output_File_Name_Path) + ".png", dpi = dpi, bbox_inches='tight')

print((datetime.now() - startTime3), ': Plot Time Taken')

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
PATH AS XYZ FORMAT OUTPUT
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Binary_Path = np.where(img2 == 1000)
Total_XYZ_1 = (Binary_Path[0]) * -1
Total_XYZ_2 = (Binary_Path[1])
        
Total_XYZ_3 = np.array([Total_XYZ_1, Total_XYZ_2])
Total_XYZ_3 = np.transpose(Total_XYZ_3)
Total_XYZ_3 = Total_XYZ_3 * cell_size
Total_XYZ_3 = np.fliplr(Total_XYZ_3)

Total_XYZ_5 = Total_XYZ_3[:,1] + max_y
Total_XYZ_5 = np.transpose(Total_XYZ_5)
Total_XYZ_4 = Total_XYZ_3[:,0] + min_x
Total_XYZ_4 = np.transpose(Total_XYZ_4)

Total_XYZ = np.array([Total_XYZ_4, Total_XYZ_5])
Total_XYZ = np.transpose(Total_XYZ)

# Creates a random z value Column, to stand out in colour
length2 = len(Total_XYZ)
ones = list(np.random.randint(1, 3, size=length2)) # Gives the z column random numbers between 1 and 3
ones = np.asarray(ones)
ones = np.reshape(ones, (length2, 1))
Total_XYZ = np.hstack((Total_XYZ, ones))

# Saving the the route as .xyz
np.savetxt("Path" + str(Output_File_Name_Path) + ".xyz", Total_XYZ, fmt='%f', delimiter=',')

# Saving the Route with its Z column in tact
empty1 = []
for n in range(len(p)):
    line = p[n]
    appended = np.where((NodalRasterWithoutBorder==line))
    empty1.append(appended)
new = np.asarray(empty1)
new = np.reshape(new, (len(new),2))
empty2 = []
for n2 in range(len(new)):
    coordsnew = new[n2]
    zvalue = WeightedRasterUntampered[coordsnew[0], coordsnew[1]]
    empty2.append(zvalue)
Total_XYZ2 = Total_XYZ
Total_XYZ2[:,2] = empty2
np.savetxt("Path with Z" + str(Output_File_Name_Path) + ".xyz", Total_XYZ2, fmt='%f', delimiter=',')
variable_pass = "{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}".format(Output_File_Name_Path, cell_size, min_x, min_y, max_x, max_y, delta_start_x, delta_start_y, delta_end_x, delta_end_y)
with open("variable pass.txt", 'w') as text_file:
      text_file.write(variable_pass)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
OUTPUT STRINGS
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
print((datetime.now() - startTime), ': TOTAL TIME')

Route_Length = length
Straight_Line_Length = (np.sqrt(abs(START_X - END_X)**2 + abs(START_Y - END_Y)**2))
print("Straight Line Length: " , str(Straight_Line_Length) + "m")
print("Route Length: " , str(Route_Length) + "m")