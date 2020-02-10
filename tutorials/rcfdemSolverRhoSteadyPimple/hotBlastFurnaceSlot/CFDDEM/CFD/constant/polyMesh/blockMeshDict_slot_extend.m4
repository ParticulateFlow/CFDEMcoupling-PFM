// Parametrized test case for a BF geometry


//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)

//Geometry

// width of wedge
m4_define(y0, 0.6)
m4_define(y1, -0.6)

// height levels
m4_define(z0, -2)
m4_define(z1, 0.75)
m4_define(z2, 3.9)
m4_define(z3, 6.4)
m4_define(z4, 24)
m4_define(z5, 28)

// xlevels
m4_define(x0, -6)
m4_define(x1, -6)
m4_define(x2, -7.4)
m4_define(x3, -7.4)
m4_define(x4, -4.7)
m4_define(x5, -4.7)


//Grid points (integers!):

m4_define(xNumberOfCells, 32)
m4_define(yNumberOfCells, 11)
m4_define(zNumberOfCells0, 25)
m4_define(zNumberOfCells1, 15)
m4_define(zNumberOfCells2, 12)
m4_define(zNumberOfCells3, 90)
m4_define(zNumberOfCells4, 10)
m4_define(rGrading, 1.0)
//m4_define(rGrading, 0.5)




/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.4.1                                 |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
(x0 y0 z0) vlabel(V0)
(0 y0 z0) vlabel(V1)
(0 y0 z1) vlabel(V2)
(x1 y0 z1) vlabel(V3)

(x2 y0 z2) vlabel(V4)
(0 y0 z2) vlabel(V5)
(0 y0 z3) vlabel(V6)
(x3 y0 z3) vlabel(V7)

(x4 y0 z4) vlabel(V8)
(0 y0 z4) vlabel(V9)
(0 y0 z5) vlabel(V10)
(x5 y0 z5) vlabel(V11)

// neg. y values
(x0 y1 z0) vlabel(V12)
(0 y1 z0) vlabel(V13)
(0 y1 z1) vlabel(V14)
(x1 y1 z1) vlabel(V15)

(x2 y1 z2) vlabel(V16)
(0 y1 z2) vlabel(V17)
(0 y1 z3) vlabel(V18)
(x3 y1 z3) vlabel(V19)

(x4 y1 z4) vlabel(V20)
(0 y1 z4) vlabel(V21)
(0 y1 z5) vlabel(V22)
(x5 y1 z5) vlabel(V23)


);

// Defining blocks:
blocks
(
    hex ( V0 V1 V2 V3 V12 V13 V14 V15 ) AB
    (xNumberOfCells zNumberOfCells0 yNumberOfCells)
    simpleGrading (rGrading 1 1)

    hex ( V3 V2 V5 V4 V15 V14 V17 V16) BC
    (xNumberOfCells zNumberOfCells1 yNumberOfCells)
    simpleGrading (rGrading 1 1)

    hex ( V4 V5 V6 V7 V16 V17 V18 V19) CD
    (xNumberOfCells zNumberOfCells2 yNumberOfCells)
    simpleGrading (rGrading 1 1)

    hex ( V7 V6 V9 V8 V19 V18 V21 V20 ) EF
    (xNumberOfCells zNumberOfCells3 yNumberOfCells)
    simpleGrading (rGrading 1 1)

    hex ( V8 V9 V10 V11 V20 V21 V22 V23 ) GH
    (xNumberOfCells zNumberOfCells4 yNumberOfCells)
    simpleGrading (rGrading 1 1)
);

// Defining patches:
boundary
(
    bottom
    {
        type patch;
        faces
        (
         (V0 V1 V13 V12)
        );
    }
    wall
    {
        type wall;
        faces
        (
         (V0 V3 V15 V12)
         (V3 V4 V16 V15)
         (V4 V7 V19 V16)
         (V7 V8 V20 V19)
         (V8 V11 V23 V20)
        );
    }
    top
    {
        type patch;
        faces
        (
         (V11 V10 V22 V23)
        );
    }
    inner
    { 
        type patch;
        faces
        (
         (V0 V1 V2 V3)
         (V3 V2 V5 V4)
         (V4 V5 V6 V7)
         (V7 V6 V9 V8)
         (V8 V9 V10 V11)

         (V12 V13 V14 V15)
         (V15 V14 V17 V16)
         (V16 V17 V18 V19)
         (V19 V18 V21 V20)
         (V20 V21 V22 V23)
        );
    }
);

mergePatchPairs 
(
);

// ************************************************************************* //
