# Oslo Model

The Oslo Model comprises of a 1-dimensional lattice of size `L`, being driven by simply adding a grain to the first site. The pile is then relaxed by determining whether each site possesses a slope above a (unique) threshold value, creating a cascading effect of grains toppling down. After each relaxation, the threshold slopes of the given site is randomly configured to either `zth = 1` or `2`.

The image below gives a visualisation of this model. Note that at the rightmost edge, grains can topple over and leave the system entirely.

![A pile of rice grains in a 1D lattice of size L. Grains can topple over the edge.](/images/rice.png?raw=true)

In order to get the results simply run the `main.py` module. This module imports 4 of the other modules and runs a script to produce the results. Make sure that all the files are kept in one folder. Remember to change the directory to the folder containing all the files so that they can be identified and read.

The `sites.py` module includes the class for the individual sites (labelled at the individual columns in the pile of grains), and `oslo.py` simulates the evolution of the pile by adding a rice grain in each iteration.

All plots and calculated values as featured in the report cannot be replicated. However, similar outputs can be produced. This is due to the random nature of the problem. It is recommended that you keep the largest power variable on `line 24` of `main.py` no larger than `7`. Also keep `N <= 10,000`, as larger values can take long to run.

NOTE: AFTER RUNNING, IT CAN TAKE UP TO 5 MINUTES FOR RESULTS TO APPEAR FOR SIZES `L <= 128` and `N <= 10,000`.

After running, the values of parameters and scaling exponents will be output on the terminal, along with 15 plots in the order of those featured in the report. Bear in mind that as system size `L=256` won't be checked, the results may not be as accurate as those in the report, where `L = 256` has been included.

## Example Results

Below is a graph showing how the height of piles for system sizes from `L = 8` to `256` evolves over time.

![Height vs time.](/images/figure_1.png?raw=true)

You can also plot the avalanche size probability against the avalanche size, `s`, itself. The data has also been binned by utilising the `log_bin_CN_2016.py` file.

![Prob vs s.](/images/figure_11.png?raw=true)

Binned data for different system sizes, `L`:

![Prob vs s for many L sizes.](/images/figure_12.png?raw=true)

These results can also be data collapsed (see `avalanchesAnalysis.py` module for the algorithm), as shown below.

![Prob vs s. data collapsed](/images/figure_16.png?raw=true)
