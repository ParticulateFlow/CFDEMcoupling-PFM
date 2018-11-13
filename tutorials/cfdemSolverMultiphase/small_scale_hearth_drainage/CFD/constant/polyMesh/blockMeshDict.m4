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
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// User-defined parameters
convertToMeters 1;

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

define(sqHalf, calc(sqSide/2))
define(sq11x, calc(-sqHalf))
define(sq11y, calc(sqHalf))
define(sq1x, calc(sqHalf))
define(sq1y, calc(sqHalf))
define(sq5x, calc(sqHalf))
define(sq5y, calc(-sqHalf))
define(sq7x, calc(-sqHalf))
define(sq7y, calc(-sqHalf))
define(c11x, calc(-cylRadius*cos(pi/4)))
define(c11y, calc(cylRadius*sin(pi/4)))
define(c1x, calc(cylRadius*cos(pi/4)))
define(c1y, calc(cylRadius*sin(pi/4)))
define(c5x, calc(cylRadius*cos(pi/4)))
define(c5y, calc(-cylRadius*sin(pi/4)))
define(c7x, calc(-cylRadius*cos(pi/4)))
define(c7y, calc(-cylRadius*sin(pi/4)))
define(ipNy, calc(cylRadius))
define(ipSy, calc(-cylRadius))
define(ipWx, calc(-cylRadius))
define(ipEx, calc(cylRadius))

vertices
(
	(sq11x sq11y 0)					//Square 11 o'clock,bottom		//0
	(sq1x sq1y 0)					//Square 1 o'clock,bottom		//1
	(sq5x sq5y 0)					//Square 5 o'clock,bottom		//2
	(sq7x sq7y 0)					//Square 7 o'clock,bottom		//3
	(c11x c11y 0)					//Cyl 11 o'clock,bottom			//4
	(c1x c1y 0)						//Cyl 1 o'clock,bottom			//5
	(c5x c5y 0)						//Cyl 5 o'clock,bottom			//6
	(c7x c7y 0)						//Cyl 7 o'clock,bottom			//7

	(sq11x sq11y cylHeight)			//Square 11 o'clock,top			//8
	(sq1x sq1y cylHeight)			//Square 1 o'clock,top			//9
	(sq5x sq5y cylHeight)			//Square 5 o'clock,top			//10
	(sq7x sq7y cylHeight)			//Square 7 o'clock,top			//11
	(c11x c11y cylHeight)			//Cyl 11 o'clock,top			//12
	(c1x c1y cylHeight)				//Cyl 1 o'clock,top				//13
	(c5x c5y cylHeight)				//Cyl 5 o'clock,top				//14
	(c7x c7y cylHeight)				//Cyl 7 o'clock,top				//15

);

blocks
(
    hex (3 2 1 0 11 10 9 8) (sqCells1D sqCells1D cellsZ) simpleGrading (1 1 1) 				//Mid-block	
    hex (7 3 0 4 15 11 8 12) (cylCellsRadial sqCells1D cellsZ) simpleGrading (1 1 1) 		//west-block 	
    hex (0 1 5 4 8 9 13 12) (sqCells1D cylCellsRadial cellsZ) simpleGrading (1 1 1) 		//north-block
    hex (2 6 5 1 10 14 13 9) (cylCellsRadial sqCells1D cellsZ) simpleGrading (1 1 1) 		//east-block
    hex (7 6 2 3 15 14 10 11) (sqCells1D cylCellsRadial cellsZ) simpleGrading (1 1 1) 		//south-block						
);

edges
(
	arc 7 4 (ipWx 0 0)					//west,bottom
	arc 4 5 (0 ipNy 0)					//north,bottom
	arc 5 6 (ipEx 0 0)					//east,bottom
	arc 6 7 (0 ipSy 0)					//south,bottom
	arc 15 12 (ipWx 0 cylHeight)		//west,top
	arc 12 13 (0 ipNy cylHeight)		//north,top
	arc 13 14 (ipEx 0 cylHeight)		//east,top
	arc 14 15 (0 ipSy cylHeight)		//south,top
);

boundary
(
    atmosphere
    {
		type patch;
		faces
        (
			(8 11 10 9)
			(8 12 15 11)
			(8 9 13 12)
			(9 10 14 13)
			(10 11 15 14)
        );
    }
    walls
    {
		type wall;
		faces
		(
		    //Bottom
			(0 1 2 3)	
			(0 3 7 4)
			(4 5 1 0)
			(2 1 5 6)
			(2 6 7 3)
			
			//West
			(15 12 4 7)
			
			//North
			(12 13 5 4)
			
			//East
			(13 14 6 5)

			//South
			(14 15 7 6)
		);
    }  
);

// ************************************************************************* //
