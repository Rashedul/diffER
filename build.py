import pybedtools
import genomepy

a = pybedtools.example_bedtool('a.bed')

windows = a.window_maker(genome ='hg38', w = 1000).saveas('bins.bed')


