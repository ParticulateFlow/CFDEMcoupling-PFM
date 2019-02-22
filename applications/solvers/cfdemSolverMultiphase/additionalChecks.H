// Additional solver-specific checks

// Useful if one wants to e.g. initialize floating particles using the Archimedes model
if (particleCloud.couplingProperties().found("unrestrictedForceModelSelection"))
{
    Warning << "Using unrestrictedForceModelSelection, results may be incorrect!" << endl;
} else
{
    #include "checkModelType.H"
}

word modelType = particleCloud.modelType();

if(!particleCloud.couplingProperties().found("useDDTvoidfraction"))
{
    Warning << "Suppressing ddt(voidfraction) is not recommended with this solver as it may generate incorrect results!" << endl;
}
