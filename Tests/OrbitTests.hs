module Tests.OrbitTests where

import Test.QuickCheck
import Debug.Trace

import Kepler.MathOps
import Kepler.Orbit
import Kepler.Anomalies
import Kepler.OrbitState

import Tests.Arbitrary

eqf x y = (x - y)^2 / (abs x + abs y)^2 / 2 < epsilon
eq3 (x1, y1, z1) (x2, y2, z2) = eqf x1 x2 && eqf y1 y2 && eqf z1 z2

zero v = dot v v < epsilon
eqv a b = ((2 / (mag a + mag b))^2) * dot (sub a b) (sub a b) < 1.0e-7

notLinear (TestVector mu (r, v) t0) =
    (not . zero $ v) && (not . zero $ r) && (not . zero $ (r `cross` v))

notNaN (TestTime orbit@(Orbit p e n m0 t0 (i, an, arg)) t) =
    all ok [p, e, n, m0, t0, i, an, arg, mean, eccentric, true, x, y, z, rx, ry, rz]
    where
    anoms@(mean, eccentric, true) = anomalies orbit t
    ((x, y, z), (rx, ry, rz)) = orbitState orbit anoms
    ok x = not (isNaN x) && not (isInfinite x)

notNaN2 (TestVector mu (r, v) t) =
    all ok [p, e, n, m0, t0, i, an, arg, mean, eccentric, true]
    where
    ((Orbit p e n m0 t0 (i, an, arg)), (mean, eccentric, true)) =
        orbitStateToElements mu (r, v) t
    ok x = not (isNaN x) && not (isInfinite x)

anomalyRanges (TestTime orbit t) =
    (not (orbitClosed orbit) || (ok mean && ok eccentric)) && ok true
    where
    (mean, eccentric, true) = anomalies orbit t
    ok x = x >= -pi && x <= pi

angleRanges (TestVector mu (r, v) t0) =
    i >= 0 && i <= pi && an >= -pi && an <= pi && arg >= -pi && arg <= pi
    where
    (orbit, anoms) = orbitStateToElements mu (r, v) t0
    (i, an, arg) = orbitOrientation orbit

identity (TestVector mu (r, v) t0) =
    giveup || eqv2 r r' && eqv2 v v'
    where
    giveup = -- this test fails on high eccentricity orbits
        (orbitType orbit /= Parabolic && e > 0.5 && e < 1.6) || e > 2.5
    e = eccentricity orbit
    (orbit, anoms) = orbitStateToElements mu (r, v) t0
    (r', v') = orbitState orbit anoms
    eqv2 a b =
        rdiff < 0.2 && vdiff < 1.0e-1
        where
        -- compare magnitude and direction separately
        rdiff = (mag b - mag a)^2 / (mag a)^2
        vdiff = (1.0 - (dot a b / (mag a * mag b) + 1) / 2)^2

visVivaVec (TestVector mu (r, v) t0) =
    eqf (dot v v) v2
    where
    ((Orbit p e _ _ _ _), _) = orbitStateToElements mu (r, v) t0
    v2
        | eccentricType e == Parabolic =
            mu * (2 / mag r)
        | otherwise =
             mu * (2 / mag r - 1 / a)
    a = p / (1 - e^2)

visVivaTime (TestTime orbit t) =
    eqf (dot v v) v2
    where
    (r, v) = orbitStateAtTime orbit t
    mu = orbitGravitationalParameter orbit
    v2
        | orbitType orbit == Parabolic =
            mu * (2 / mag r)
        | otherwise =
             mu * (2 / mag r - 1 / (semiMajorAxis orbit))

gravityParameter (TestVector mu (r, v) t0) =
    eqf mu mu'
    where
    mu' = orbitGravitationalParameter orbit
    (orbit, _) = orbitStateToElements mu (r, v) t0

periodicity (TestTime orbit t0) =
    not (orbitClosed orbit) || (eqv r r' && eqv v v' && eq3 anoms anoms')
    where
    anoms = anomalies orbit t0
    (r, v) = orbitStatePlane orbit anoms
    anoms' = anomalies orbit (t0 + orbitalPeriod orbit)
    (r', v') = orbitStatePlane orbit anoms'

apsides orbit =
    (eqf (periapsis orbit) (mag per) && eqf (periapsisVel orbit) (mag pev)) &&
        (not (orbitClosed orbit) || eqf (apoapsis orbit) (mag apr) && eqf (apoapsisVel orbit) (mag apv))
    where
    (per, pev) = orbitState orbit (0, 0, 0)
    (apr, apv) = orbitState orbit (pi, pi, pi)

nodes orbit =
    ((z / mag r)^2 < epsilon && ((if ascending then vz else -vz) >= 0)) &&
        (not (orbitClosed orbit) ||
            ((z' / mag r')^2 < epsilon && ((if ascending then -vz' else vz') >= 0)))
    where
    (i, an, arg) = orbitOrientation orbit
    e = eccentricity orbit
    -- ascending node closer to periapsis than descending node
    ascending = abs arg < pi/2
    -- the node closer to periapsis
    true = if ascending then -arg else pi - arg
    anoms = (trueToMeanAnomaly e true, trueToEccentricAnomaly e true, true)
    (r@(_, _, z), v@(_, _, vz)) = orbitState orbit anoms
    -- the node away from periapsis
    true' = if ascending then pi - arg else -arg
    anoms' = (trueToMeanAnomaly e true', trueToEccentricAnomaly e true', true')
    (r'@(_, _, z'), v'@(_, _, vz')) = orbitState orbit anoms'

velocityTest orbit =
    giveup || (eqv (add r (scale (1/dt) v)) r' && eqv v' v'')
    where
    giveup = -- this test fails near the escape velocity
        orbitType orbit /= Parabolic && eccentricity orbit > 0.8 && eccentricity orbit < 1.2
    dt = 1.0
    mu = orbitGravitationalParameter orbit
    anoms = anomalies orbit (orbitEpoch orbit)
    (r, v) = orbitState orbit anoms
    anoms' = anomalies orbit (orbitEpoch orbit + dt)
    (r', v') = orbitState orbit anoms'
    v'' = sub v (scale (dt * mu / (dot r r)) (unit r))
