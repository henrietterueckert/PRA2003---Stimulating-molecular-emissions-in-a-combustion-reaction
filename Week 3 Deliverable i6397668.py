# IMPORTANT: Make sure you have any of the output files in your directory
# This code calculates the average number of a chosen particle type per event (500,000 runs) with the uncertainty

# Equations i used:
# N = number of events
# mean = (sum of counts) / N
# variance = sum((x - mean)^2) / (N - 1)
# uncertainty    = sqrt(variance / N)   

# important for this code: Python always counts from 0, so e.g. if ID is in column 4, then i have to put in column 3

# Code start
import math  # needed for sqrt() when computing the uncertainty

# Ask the user for the inputs instead of hardcoding them
filename = input("Enter the data file name (like output-Setx.txt): ") # same script works for any data file

# Print ID choosing options
print('Possible particle IDs to chooose:  211, -211, 321, -321, 2212, -2212, 3122, -3122, 3312, -3312, 3334, -3334')

# Crash protection
try:
    id_of_interest = int(input("Enter the particle ID to count: ")) # same script works for any particle ID
except ValueError:  # int() raises ValueError if the user doesn't type a number
    print("Error: the ID must be an integer")
    exit()

# Check the file exists before using it. Protection
try:
    f = open(filename, "r")
except FileNotFoundError:
    print("Error: could not find the file", filename)
    exit()

#  keep one count per event in a list
counts = []   # counts[i] = n. of matching particles found in event i

# loop reads 1 event with header line and the particle lines that belong to it
while True:
    header_line = f.readline()
    if header_line == "":   # means there are no lines left to read!
        break

    header = header_line.split()   # split the line between spaces

    # A header should have 2 things onnly!! (event number and particle count)
    if len(header) != 2: # if more than 2 it will stop reading
        print("Error: malformed header line:", header_line)
        break

    # Particle counts
    try:
        n_particles = int(header[1]) # the number of particles converted to integer, tells the loop how many lines in event
    except ValueError:
        print("Error, count is not an integer", header_line)
        break

# Initialising
    count = 0   # n. particle of interest in an event

    for i in range(n_particles): # reads n_particle lines, next readline will be next event
        line = f.readline()

# Making sure each particle line has 4 components
        columns = line.split() # spltting components between spaces
      
        if len(columns) != 4: 
            print("Error: malformed data line:", line)
            continue # doesn't break, it will skip this line!!

        try:
            particle_id = int(columns[3]) # convert ID to integer
        except ValueError:
            print("Error: ID is not an integer:", line)
            continue # again if this fails it will only skip a line instead of breaking

        # counting particle of intest
        if particle_id == id_of_interest: # if clause to only tally the particl ID that matches the one we chose
            count = count + 1 # tallies up

# Loop is finished for one event. add this event's count to the list
    counts.append(count)

f.close()   # closing file

# Number of events read
n_events = len(counts)

# Average number of matching particles per event
if n_events == 0: # cannot calculate mean from 0
    print("division by 0 error.")
    exit()

mean = sum(counts) / n_events


# The sample variance divides by (N - 1) so n_events needs to be >1
if n_events > 1:
    variance = sum((x - mean) ** 2 for x in counts) / (n_events - 1)   # sample variance equation
    error = math.sqrt(variance / n_events)   # uncertainty equation
else:
    error = 0.0 # there is no variance otherwise

# Printing result
print("Events read:", n_events)
print("Average per event:", round(mean, 4), "+/-", round(error, 4)) # round() to 4 decimal places keeps the output readable