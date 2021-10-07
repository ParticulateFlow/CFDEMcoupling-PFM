# Test large, incomplete data bases

This tutorial tests the use of data bases (in memory), which are smaller
than the number of snaphots in the data base (on disk). This distinction
becomes important for cases with large meshes, a large number of snapshots,
or cases with both of these properties. In such cases, the memory (RAM) of
the computer might be too small to fit all snapshots.

In such a case - with *M*, the number of snapshots in memory, being smaller
than *N*, the number of snapshots on disk - we use the recurrence model
designed for large, complete data bases.

This tutorial is based on the *laminar flow over cylinder* tutorial. It is
rather small, but is nevertheless used to demonstrate the use of the large,
incomplete data bases feature.


## Case description

This case is run for a small amount of time, and then, the tool `rStatAnalysis`
is used as a post-processing tool to compute the recurrence matrix  [1].
In this case, `rStatAnalysis` is run with a data base (in memory) smaller than
the number of snapshots on disk.


## Results

If we monitor the log file written by `rStatAnalysis`, we fill find the line
similar to the following:

`Checking all 20 fields`

This line tells us the total number of snapshots on disk, which will all be
checked for existence and readabiltiy.

Afterwards, when it is enured that all expected fields are present, the data
base in memory is filled with fields read from disk. However, as the data base
in memory is smaller than the number of snapshots on disk, the data base in
memory will be incomplete.

Thus, we will encounter a line similar to the following one, which tells us
that a number of fields is being read.

`Reading fields for 11 dataBase slots`

When computing the recurrence matrix, all *N* snapshots need to be compared
with each other. As the data base holds only *M* snapshots, with *M<N*, at
some point, the snapshots in the data base will need to be replaced with
other ones, wich have not been read at that point.
Such an operation is indicated by a line as shown below, provided the switch
`verboseVacate` is enabled in `recProperties`.

` --> vacate dataBase  : 5`

After the recurrence matrix has been computed, the recurrence model prints
the total number of reads-from-disk to the Terminal. This is shown below:

`Nr. of reads from disk : 30; compared to a theoretical minimum of 20 reads.`

This indicates, that the total number of reads-from-disk is necessarily
larger that the number of snapshots. The total number of reads-from-disk
depends on the relative sizes of *N* and *M*. However, it will always be larger
than the theoretical minimum, which is achieved, when *M* is equal to *N*, i.e.
each field is read only once.


## Tested

This collection of cases has been tested with:

* OpenFOAM-5.0


## References

[1] T. Lichtenegger and S. Pirker. Recurrence CFD – A novel approach to simulate
multiphase flows with strongly separated time scales. Chemical Engineering Science,
153:394-410, 2016.
