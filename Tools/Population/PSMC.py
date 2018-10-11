import os
import numpy as np
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()

from Tools.Abstract import Tool


class PSMC(Tool):
    def __init__(self, path="", max_threads=4, max_memory="100G", max_per_thread_memory="5G"):
        Tool.__init__(self, "plink", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

    def psmc_plot(self, sample_label_list, psmc_list, generation_time, absolute_mutation_rate,
                  output_prefix, plot_grid=False):
        """
        -u FLOAT   absolute mutation rate per nucleotide [2.5e-08]
         -s INT     skip used in data preparation [100]
         -X FLOAT   maximum generations, 0 for auto [0]
         -x FLOAT   minimum generations, 0 for auto [10000]
         -Y FLOAT   maximum popsize, 0 for auto [0]
         -m INT     minimum number of iteration [5]
         -n INT     take n-th iteration (suppress GOF) [20]
         -M titles  multiline mode [null]
         -f STR     font for title, labels and tics [Helvetica,16]
         -g INT     number of years per generation [25]
         -w INT     line width [4]
         -P STR     position of the keys [right top]
         -T STR     figure title [null]
         -N FLOAT   false negative rate [0]
         -S         no scaling
         -L         show the last bin
         -p         convert to PDF (with epstopdf)
         -R         do not remove temporary files
         -G         plot grid
        """

        options = " -g %i" % generation_time
        options += " -u %f" % absolute_mutation_rate
        options += " -p" # TODO: find sence of this options
        options += " -G" if plot_grid else plot_grid
        options += " -M %s" % ",".join(sample_label_list)
        options += " %s" % output_prefix
        options += " %s" % " ".join(psmc_list)

        self.execute(options=options, cmd="psmc_plot.pl")
