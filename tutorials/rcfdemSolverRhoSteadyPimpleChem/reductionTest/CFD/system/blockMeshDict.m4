/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           | Copyright (C) 2016 Ehsan Madadi-Kandjani        |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General macros to create cylinder mesh

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(btQuad, ($1b $2b $2t $1t))
define(topQuad, ($1t $4t $3t $2t))
define(bottomQuad, ($1b $2b $3b $4b))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1.0;

// Inner square side half
define(s, 0.02)

// Inner square side curvature
define(sc, 0.025)

// cylinder radius
define(r, 0.05)

// Height of cylinder
define(z, 0.2)

// Base z
define(Zb, 0)

// Outlet z
define(Zt, calc(Zb + z))

// Number of cells at inner square
define(Ns, 6)

// Number of cells between inner square and circle
define(Ni, 4)

// Number of cells in the cylinder height
define(Nz, 20)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

define(vert, (x$1$2 y$1$2 $3))
define(evert, (ex$1$2 ey$1$2 $3))

// 45 degree points angle
define(a0, -45)
define(a1, -135)
define(a2, 135)
define(a3, 45)

// Half of 45 degree points angle
define(ea0, 0)
define(ea1, -90)
define(ea2, 180)
define(ea3, 90)

define(ca0, calc(cos((pi/180)*a0)))
define(ca1, calc(cos((pi/180)*a1)))
define(ca2, calc(cos((pi/180)*a2)))
define(ca3, calc(cos((pi/180)*a3)))

define(sa0, calc(sin((pi/180)*a0)))
define(sa1, calc(sin((pi/180)*a1)))
define(sa2, calc(sin((pi/180)*a2)))
define(sa3, calc(sin((pi/180)*a3)))

define(cea0, calc(cos((pi/180)*ea0)))
define(cea1, calc(cos((pi/180)*ea1)))
define(cea2, calc(cos((pi/180)*ea2)))
define(cea3, calc(cos((pi/180)*ea3)))

define(sea0, calc(sin((pi/180)*ea0)))
define(sea1, calc(sin((pi/180)*ea1)))
define(sea2, calc(sin((pi/180)*ea2)))
define(sea3, calc(sin((pi/180)*ea3)))

// Inner square x and y position

// x
define(x00, s)
define(x01, calc(-1.0*s))
define(x02, calc(-1.0*s))
define(x03, s)

// y
define(y00, calc(-1.0*s))
define(y01, calc(-1.0*s))
define(y02, s)
define(y03, s)

// Circle x and y positions

// x
define(x10, calc(r*ca0))
define(x11, calc(r*ca1))
define(x12, calc(r*ca2))
define(x13, calc(r*ca3))

// y
define(y10, calc(r*sa0))
define(y11, calc(r*sa1))
define(y12, calc(r*sa2))
define(y13, calc(r*sa3))

// Inner square x and y position middle curvatures

// x
define(ex00, sc)
define(ex01, 0)
define(ex02, calc(-1.0*sc))
define(ex03, 0)

// y
define(ey00, 0)
define(ey01, calc(-1.0*sc))
define(ey02, 0)
define(ey03, sc)

// Circle x and y positions middle curvatures

// x
define(ex10, calc(r*cea0))
define(ex11, calc(r*cea1))
define(ex12, calc(r*cea2))
define(ex13, calc(r*cea3))

// y
define(ey10, calc(r*sea0))
define(ey11, calc(r*sea1))
define(ey12, calc(r*sea2))
define(ey13, calc(r*sea3))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(s0b)
    vert(0, 1, Zb) vlabel(s1b)
    vert(0, 2, Zb) vlabel(s2b)
    vert(0, 3, Zb) vlabel(s3b)
    
    vert(1, 0, Zb) vlabel(r0b)
    vert(1, 1, Zb) vlabel(r1b)
    vert(1, 2, Zb) vlabel(r2b)
    vert(1, 3, Zb) vlabel(r3b)
    
    vert(0, 0, Zt) vlabel(s0t)
    vert(0, 1, Zt) vlabel(s1t)
    vert(0, 2, Zt) vlabel(s2t)
    vert(0, 3, Zt) vlabel(s3t)
    
    vert(1, 0, Zt) vlabel(r0t)
    vert(1, 1, Zt) vlabel(r1t)
    vert(1, 2, Zt) vlabel(r2t)
    vert(1, 3, Zt) vlabel(r3t)
);

blocks
(
    //block0
    hex2D(s1, s0, s3, s2)
    square
    (Ns Ns Nz)
    simpleGrading (1 1 1)
    
    //block1
    hex2D(s0, r0, r3, s3)
    innerCircle
    (Ni Ns Nz)
    simpleGrading (0.5 1 1)
    
    //block2
    hex2D(s3, r3, r2, s2)
    innerCircle
    (Ni Ns Nz)
    simpleGrading (0.5 1 1)
    
    //block3
    hex2D(s2, r2, r1, s1)
    innerCircle
    (Ni Ns Nz)
    simpleGrading (0.5 1 1)
    
    //block4
    hex2D(s1, r1, r0, s0)
    innerCircle
    (Ni Ns Nz)
    simpleGrading (0.5 1 1)
);

edges
(
    //Circle edges
    arc r3b r0b evert(1, 0, Zb)
    arc r0b r1b evert(1, 1, Zb)
    arc r1b r2b evert(1, 2, Zb)
    arc r2b r3b evert(1, 3, Zb)
    
    //Circle edges
    arc r3t r0t evert(1, 0, Zt)
    arc r0t r1t evert(1, 1, Zt)
    arc r1t r2t evert(1, 2, Zt)
    arc r2t r3t evert(1, 3, Zt)
    
    arc s3b s0b evert(0, 0, Zb)
    arc s0b s1b evert(0, 1, Zb)
    arc s1b s2b evert(0, 2, Zb)
    arc s2b s3b evert(0, 3, Zb)
    
    arc s3t s0t evert(0, 0, Zt)
    arc s0t s1t evert(0, 1, Zt)
    arc s1t s2t evert(0, 2, Zt)
    arc s2t s3t evert(0, 3, Zt)
    
);

patches
(
    wall walls
    (
        btQuad(r0, r3)
        btQuad(r1, r0)
        btQuad(r2, r1)
        btQuad(r3, r2)
    )

    patch inlet
    (
        bottomQuad(s3, s0, s1, s2)
        bottomQuad(s3, r3, r0, s0)
        bottomQuad(s2, r2, r3, s3)
        bottomQuad(s1, r1, r2, s2)
        bottomQuad(s0, r0, r1, s1)
    )
    
    patch outlet
    (
        topQuad(s3, s0, s1, s2)
        topQuad(s3, r3, r0, s0)
        topQuad(s2, r2, r3, s3)
        topQuad(s1, r1, r2, s2)
        topQuad(s0, r0, r1, s1)
    )
);

mergePatchPairs
(
);
