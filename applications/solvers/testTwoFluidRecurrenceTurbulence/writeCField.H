/* ----------------------------------------------
	Write averaged particle volume fraction
---------------------------------------------- */

// essential information
Info << "Total mass :";
Info << tab << sum(C*rhoCarrier*alphaCarrier*mesh.V());
Info << endl;

Info << "Total Carrier mass :";
Info << tab << sum(rhoCarrier*alphaCarrier*mesh.V());
Info << endl;
