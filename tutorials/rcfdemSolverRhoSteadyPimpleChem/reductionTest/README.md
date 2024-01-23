## Test case for ore reduction in a small cylinder

1. fill a cylinder with sperical particles (initCase.sh)
2. create CFD mesh (mesh.sh)
3. introduce gas flow to reduce ore (runCase.sh)


#### Comments:

After particles have settled, they remain frozen to their positions. The mass loss over time is recorded in DEM/monitor/mass_layer.dat [in kg] and can be compared to the target value target_mass_loss_pellets.txt [in g].
