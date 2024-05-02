# build genome database

import pybedtools
pybedtools.bedtool.BedTool.window_maker(...)

from pybedtools import BedTool

# Example bed files
bed_file1 = "file1.bed"
bed_file2 = "file2.bed"

# Creating BedTool objects from bed files
bed1 = BedTool(bed_file1)
bed2 = BedTool(bed_file2)

# Defining the window size
window_size = 1000

# Making windows around bed1
windows = bed1.window_maker(g=bed2, w=window_size)

# Iterating through windows
for window in windows:
    print(window)
