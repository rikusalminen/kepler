module Kepler.MathOps (
    epsilon,
    inf,
    add,
    sub,
    mul,
    scale,
    dot,
    cross,
    mag,
    unit,
    vproduct,
    sign,
    clamp,
    angleClamp
    ) where


epsilon = 1.0e-10 :: Double
inf = 1.0/0.0 :: Double

add (x1, y1, z1) (x2, y2, z2) = (x1+x2, y1+y2, z1+z2)
sub (x1, y1, z1) (x2, y2, z2) = (x1-x2, y1-y2, z1-z2)
mul (x1, y1, z1) (x2, y2, z2) = (x1*x2, y1*y2, z1*z2)
scale c (x1, y1, z1) = (c*x1, c*y1, c*z1)

dot (x1, y1, z1) (x2, y2, z2) = x1*x2 + y1*y2 + z1*z2
cross (x1, y1, z1) (x2, y2, z2) = (y1*z2 - z1*y2, -(x1*z2 - z1*x2), x1*y2 - y1*x2)

mag v = sqrt (dot v v)
unit v = scale (1/mag v) v

vproduct (x, y, z) v = (dot x v, dot y v, dot z v)

sign x = if x < 0 then -1 else 1
clamp mi ma x = min ma (max x mi)

angleClamp x
    | abs x <= pi =
        x
    | otherwise =
        -pi+2*pi*fract((x+pi)/(2*pi))
    where
    fract f = f - (fromIntegral $ floor f)
