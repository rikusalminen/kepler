module Kepler.OrbitIntercept (
    OrbitIntercept(..),
    orbitIntercept
    ) where

import Kepler.MathOps
import Kepler.Orbit
import Kepler.Anomalies
import Kepler.OrbitState

import Debug.Trace

data OrbitIntercept = OrbitIntercept {
    interceptTime :: Double,
    interceptDistance :: (Double, Double, Double),
    interceptVelocity :: (Double, Double, Double)
    }
    deriving Show

orbitIntercept :: (Orbit, Orbit) -> (Double, Double) -> Double -> [OrbitIntercept]
orbitIntercept (o1, o2) (t0, t1) distance
    | (mu1-mu2)^2/(mu1*mu2) > epsilon =
        error "Orbits around different bodies don't intercept"
    | t0 >= t1 =
        []
    | orbitClosed o1 && orbitalPeriod o1 / 2 - (t1 - t0) < -epsilon =
        orbitIntercept (o1, o2) (t0, t0 + orbitalPeriod o1 / 2) distance ++
        orbitIntercept (o1, o2) (t0 + orbitalPeriod o1 / 2, t1) distance
    | orbitClosed o2 && orbitalPeriod o2 / 2 - (t1 - t0) < -epsilon =
        orbitIntercept (o1, o2) (t0, t0 + orbitalPeriod o2 / 2 ) distance ++
        orbitIntercept (o1, o2) (t0 + orbitalPeriod o2 / 2, t1) distance
    | m1t0 <= 0 && m1t1 >= 0 && m1t0^2 > epsilon && m1t1^2 > epsilon && (m1t1 - m1t0)^2 > epsilon =
        let t = t0 - m1t0 * (t1 - t0) / (m1t1 - m1t0)
        in
            orbitIntercept (o1, o2) (t0, t) distance ++
            orbitIntercept (o1, o2) (t, t1) distance

    | m2t0 <= 0 && m2t1 >= 0 && m2t0^2 > epsilon && m2t1^2 > epsilon && (m2t1 - m2t0)^2 > epsilon =
        let t = t0 - m2t0 * (t1 - t0) / (m2t1 - m2t0)
        in
            orbitIntercept (o1, o2) (t0, t) distance ++
            orbitIntercept (o1, o2) (t, t1) distance
    | otherwise =
        orbitIntercept' (o1, o2) (t0, t1) distance maxSplits $
            ((orbitStateAtTime o1 t0, orbitStateAtTime o2 t0), (orbitStateAtTime o1 t1, orbitStateAtTime o2 t1))
    where
    (mu1, mu2) = (orbitGravitationalParameter o1, orbitGravitationalParameter o2)
    (m1t0, m1t1) = (meanAnomalyAtTime o1 t0, meanAnomalyAtTime o1 t1)
    (m2t0, m2t1) = (meanAnomalyAtTime o2 t0, meanAnomalyAtTime o2 t1)
    maxSplits = 1

-- Finds at most 2^splits interceptions
orbitIntercept' (o1, o2) (t0, t1) distance splits (s0@((r1t0, v1t0), (r2t0, v2t0)), s1@((r1t1, v1t1), (r2t1, v2t1)))
    | rv0^2 < epsilon || (dot d0 d0) < epsilon =
        traceShow ("hit begin", t0) $
        [OrbitIntercept t0 d0 v0]
        --orbitIntercept' (o1, o2) (t0 + timestep, t1) distance splits
    | distance > 0 && (mag d0 - distance)^2 < epsilon =
        traceShow ("hit begin", t0) $
        [OrbitIntercept t0 d0 v0]
    | t1 - t0 < timestep && distance > 0.0 =
        if (mag d0 - distance) * (mag d1 - distance) <= 0
            then [OrbitIntercept tm dm vm]
            else []
    | t1 - t0 < timestep =
        if rv0 * rv1 <= 0.0
        --if rv0 >= 0.0 && rv1 <= 0.0
            then [OrbitIntercept tm dm vm]
            else []
    | distance > 0.0 && (mag d0 - distance) * (mag d1 - distance) > 0 =
        -- split
        if splits > 0
            then left (splits-1) ++ right (splits-1)
            else []
    | distance > 0.0 =
        -- choose
        if (mag d0 - distance) * (mag dm - distance) <= 0 then left splits else right splits
    | rv0 * rv1 > 0 =
        -- split
        trace "split" $
        if splits > 0
            then left (splits-1) ++ right (splits-1)
            else []
    | otherwise =
        -- choose
        traceShow (r1t0, r2t0) $
        traceShow (r1t1, r2t1) $
        traceShow (v1t0, v2t0) $
        traceShow (v1t1, v2t1) $
        traceShow (mag d0, mag v0) $
        traceShow (rv0, rv1) $
        trace "choose" $
        if rv0 * rvm <= 0.0 then left splits else right splits
    where
    (d0, d1) = (sub r2t0 r1t0, sub r2t1 r1t1)
    (v0, v1) = (sub v1t0 v2t0, sub v1t1 v2t1)
    (rv0, rv1) = (v0 `dot` unit d0, v1 `dot` unit d1)
    tm = (t0 + t1) / 2 -- TODO: better heuristic than binary search!
    sm@((r1tm, v1tm), (r2tm, v2tm)) = (orbitStateAtTime o1 tm, orbitStateAtTime o2 tm)
    (dm, vm) = (sub r2tm r1tm, sub v1tm v2tm)
    rvm = vm `dot` unit dm
    left s = orbitIntercept' (o1, o2) (t0, tm) distance s (s0, sm)
    right s = orbitIntercept' (o1, o2) (tm, t1) distance s (sm, s1)
    timestep = 1.0

hillclimb (r1, v1, a1) (r2, v2, a2) =
    -- newton-rhapson step
    -r / dr
    where
    -- r = d/dt (r1-r2)^2
    r = (mag r1 * mag v1) + (mag r2 * mag v2) - (dot v1 r2) - (dot r1 v2)
    -- dr = d/dt r
    dr =
        (dot v1 v1 + mag r1 * mag a1) +
        (dot v2 v2 + mag r2 * mag a2) -
        (dot a1 r2) -
        (dot r1 a2) -
        2 * (dot v1 v2)

{-
orbitIntercept (o1, o2) (t0, t1)
    | t0 > t1 =
        []
    | (rv0 < 0 || rv1 < 0) || (rv0 * rv1 > 0) =
        if (t1 - t0) <= min_interval
            then []
            else
                orbitIntercept (o1, o2) (t0, t0+t1/2) ++
                orbitIntercept (o1, o2) (t0+t1/2, t1)
    | otherwise =
        []
    where
    mu = orbitGravitationalParameter o1
    gravity r = (-1.0 * mu / (mag r)^3) `scale` r

    ((r1t0, v1t0), (r2t0, v2t0)) = (orbitStateAtTime o1 t0, orbitStateAtTime o2 t0)
    ((r1t1, v1t1), (r2t1, v2t1)) = (orbitStateAtTime o1 t1, orbitStateAtTime o2 t1)
    (a1t0, a2t0) = (gravity r1t0, gravity r2t0)
    (a1t1, a2t1) = (gravity r1t1, gravity r2t1)

    (d0, d1) = (sub r2t0 r1t0, sub r2t1 r1t1)
    (v0, v1) = (sub v1t0 v2t0, sub v1t1 v2t1)
    (rv0, rv1) = (v0 `dot` unit d0, v1 `dot` unit d1)

    t0' = t0 + (hillclimb (r1t0, v1t0, a1t0) (t2t0, v2t0, a2t0))
    t1' = t1 + (hillclimb (r1t1, v1t1, a1t1) (t2t1, v2t1, a2t1))

    dt = 1.0
    min_interval = 60.0
    -}

testOrbit mu ap pe =
    Orbit p e n 0 0 (0, 0, 0)
    where
    a = (ap + pe) / 2.0
    e = (ap - pe) / (ap + pe)
    p = a * (1 - e^2)
    n = sqrt (mu / a^3)

test =
    traceShow ("t", t) $
    traceShow ("dr", dr) $
    traceShow ("dv", dv) $
    traceShow ("mut", orbitGravitationalParameter ot) $
    traceShow ("mu2", orbitGravitationalParameter o2) $
    traceShow ("ot", ot) $
    traceShow ("o2", o2) $
    traceShow ("ot", periapsis ot, apoapsis ot, eccentricity ot, orbitalPeriod ot) $
    traceShow ("o2", periapsis o2, apoapsis o2, eccentricity o2, orbitalPeriod o2) $
    traceShow ("ot", orbitStateAtTime ot t) $
    traceShow ("o2", orbitStateAtTime o2 t) $
    interceptions
    where
    -- textbook hohmann transfer orbit
    mu = 40.0e13
    o1 = (testOrbit mu 7000.0e3 7000.0e3)
    o2 = (testOrbit mu 7500.0e3 7500.0e3) { meanAnomalyAtEpoch = m }
    ot = (testOrbit mu 7560.0e3 7000.0e3)
    --pt = orbitalPeriod ot / 2.0
    pt = 3066
    m = pi - pt * meanMotion o2
    interceptions = orbitIntercept (ot, o2) (-1.0, 1.0 * orbitalPeriod o2) (-1.0)
    (OrbitIntercept t dr dv) = head interceptions

main = do
    print [(interceptTime i, mag . interceptDistance $ i, mag . interceptVelocity $ i) | i <- test]

