# Feedback

## Goal of the project
The goal of the project seems clear (simulating F1 hybrid genomes) and has an appropriate plan. It may help to automatically include pandas and matplotlib as dependencies to simplify the decision-making for the user. But otherwise, this is straightforward and well-thought out.

## The data
The required input data files are well-explained, and the resulting output options are also given. It reads as though there will always be three output files (VCF, CSV, and plots) and one optional stats file. It may be useful to give the user the option of specifying which output files they want (i.e., one is necessary, but all are optional) for more customization.

## The code
- The current code/pseudocode in `distortopia/segdistorters.py` and `distortopia/simf1poly.py` could use more comments do lay out what the functions are doing and how they fit into the overall picture. Making a more detailed skeleton (even just comments) will help guide the project.
- Currently, the code is in two separate python script files, with a different function in each file. Both are imported into the `__main__.py` file for use. If it is run (I think) either one of the functions is executed to simulate F1s or detect distortion.
- Individual functions that can help accoplish part of the project would be:
	- Write a CSV of observed vs. expected frequencies (could be broken into observed function and expected function, and then a small function to combine them and write the CSV).
	- Write a function to plot the distorted regions, beginning with one default plot. Later, it can be built out with different options for different types of plots.


## Code contributions/ideas
In the `segdistorters.py` file, I added a toggle for writing the observed and expected ratios to a CSV. I haven't tested it, so I don't know if it shows up as just lists, or if the values are properly broken into different rows. If you change the default to `None` rather than `False`, and change `if csv` to `if type(csv) is str`, then you can have the user specify the name of the resulting CSV file instead of the default name I gave it.