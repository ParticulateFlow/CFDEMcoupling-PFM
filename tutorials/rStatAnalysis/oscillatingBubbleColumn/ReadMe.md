# oscillatingBubbleColumn

## Case description

This case demonstrates the use of the utility `rStatAnalysis` to create the recurrence 
statistics and create a plot of the recurrence matrix. The underlying simulation case 
is a plain, quasi 2D bubble column, which is computed by `reactingTwoPhaseEulerFoam`.

After the two-fluid solver has finished, `rStatAnalysis` is run to create the 
recurrence statistics. A slight modification to the class `standardRecModel`, now 
allows to specify an alternative path to the recurrence data base, which defaults 
to a sub-directory of the case-directory named *dataBase*. In this case, however, 
we specify the case-directory to act as the data base, as be use `rStatAnalysis` 
purely as a post-processing tool for the current case.

## Tested

This case has been tested with:

* OpenFOAM-4.1
