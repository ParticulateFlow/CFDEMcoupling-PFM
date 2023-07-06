## Blast furnace slot model including heat transfer

1. First, the BF geometry needs to be filled and the layer structure needs to be created.
   To this end, execute run_init.sh. The resulting data is transferred to the appropriate places by post_init.sh which is automatically called by run_init.sh.

2. Then, a database with the mean flow fields needs to be created. First, execute run_CFDDEM.sh.

3. Once finished, a data-driven run can be performed. Execute run_dataDrivenCFD.sh.

Repeat steps (2) and (3) as often as desired.
