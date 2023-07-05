Test case for sticky fines transport through packed particle column

1) fill a cylinder with sperical particles (parDEMrun.sh)
2) create CFD mesh (mesh.sh)
3) introduce gas flow with fines (parCFDDEMrun.sh)


Comments:

After particles have settled, they are removed from the bottom region so that a moving bed is created.
Fines are deposited because of their stickiness if the local flow velocity is below a critical value to be set in CFD/dustParameters.
