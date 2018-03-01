/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// User-defined parameters
define(cylRadius, 0.1)
define(cylHeight, 0.08)
define(sqSide, 0.1)
define(outletSize, 0.02)
define(outletStartZ, 0.04)
define(outletAngularPos, 2.5) //degrees from x-axis
define(sqCells1D, 8)
define(cylCellsRadial, 3)
define(cellsZ, 6)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Derived parameters
define(midPointX, calc(cylRadius*cos(outletAngularPos*pi/180)))
define(midPointY, calc(cylRadius*sin(outletAngularPos*pi/180)))

actions
(
	{
		name    c_out;
		type    cellSet;
		action  new;
		source  boxToCell;
		sourceInfo
		{
		   box (calc(midPointX-0.5*outletSize) calc(midPointY-0.5*outletSize) outletStartZ)(calc(midPointX+0.5*outletSize) calc(midPointY+0.5*outletSize) calc(outletStartZ+1.0*outletSize));
		}
 	}

 	{
		name    outlet;
		type    faceSet;
		action  new;
		source  patchToFace;
		sourceInfo
		{
		   name "wall";
		}
  	}

 	{
		name    outlet;
		type    faceSet;
		action  subset;
		source  cellToFace;
		sourceInfo
		{
		   set c_out;
		   option all;
		}
 	}
);
// ************************************************************************* //
