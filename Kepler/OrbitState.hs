module Kepler.OrbitState (
    orbitStateAtTime,
    orbitState,
    orbitStatePlane,
    orbitStateToElements
    ) where

import Kepler.MathOps (epsilon, sub, scale, dot, cross, mag, unit, vproduct, sign, clamp)
import Kepler.Orbit (Orbit(..), OrbitType(..), orbitMatrix, eccentricType)
import Kepler.Anomalies (trueToEccentricAnomaly, eccentricToMeanAnomaly, anomalies)

orbitStateAtTime orbit = orbitState orbit . anomalies orbit

orbitState orbit anoms =
    (vproduct matrix r, vproduct matrix v)
    where
    matrix = orbitMatrix orbit
    (r, v) = orbitStatePlane orbit anoms

orbitStatePlane (Orbit p e n _ _ _) (mean, eccentric, true) =
    state (eccentricType e)
    where
    state Circular =
        ((x, y, 0), (vx, vy, 0))
        where
        x = p * cos eccentric
        y = p * sin eccentric
        vx = -v * sin eccentric
        vy = v * cos eccentric
        v = p * n
    state Parabolic =
        ((x, y, 0), (vx, vy, 0))
        where
        r = p / (1 + cos true)
        x = r * cos true
        y = r * sin true
        dx = (-p * sin true) / (1 + cos true)^2
        dy = p / (1 + cos true)
        d = sqrt (dx^2 + dy^2)
        vx = v * dx / d
        vy = v * dy / d
        v = n * sqrt (p^2 / 2 * (1 + cos true))
    state Elliptic =
        ((x, y, 0), (vx, vy, 0))
        where
        b = a * sqrt (1 - e^2)
        x = a * (cos eccentric - e)
        y = b * sin eccentric
        edot = n / (1 - e * cos eccentric)
        vx = -a * sin eccentric * edot
        vy =  b * cos eccentric * edot
    state Hyperbolic =
        ((x, y, 0), (vx, vy, 0))
        where
        --x = -a * (e - cosh eccentric)
        --y = -a * sqrt(e^2 - 1) * sinh eccentric
        r = p / (1 + e * cos true)
        x = r * cos true
        y = r * sin true
        -- mu = (n^2 * (p / abs (1 - e^2))^3)
        mu = (n^2 * (-p / (1 - e^2))^3)
        --v = sqrt ((mu / p) * (1 + e^2 + 2 * cos true))
        v = sqrt (mu * ((2 / r) - (1 / a)))
        dx = -p * sin true / (1 + e * cos true)^2
        dy = p * (e + cos true) / (1 + e * cos true)^2
        d = sqrt (dx^2 + dy^2)
        vx = v * dx / d
        vy = v * dy / d
    a = p / (1 - e^2)

orbitStateToElements mu (r, v) t0 =
    trajectory
    where
    trajectory
        | mag v < epsilon || mag r < epsilon || mag k < epsilon =
            linearTrajectory
        | otherwise =
            conicTrajectory (eccentricType e)
    linearTrajectory =
        undefined

    -- specific relative angular momentum, direction = orbital plane normal
    k = cross r v
    -- eccentricity vector, length = eccentricity, direction = to periapsis
    evec = sub (scale (-1/mu) (cross k v)) (unit r)
    e = mag evec
    -- line of nodes, length = zero for equatorial, direction = to ascending node
    -- nodes = cross (0, 0, 1) k
    nodes = let (kx, ky, _) = k in (-ky, kx, 0)

    conicTrajectory Circular =
        (Orbit radius e n true t0 (i, an, arg), (true, true, true))
        where
        radius = mu / (dot v v)
        n = sqrt (mu / radius^3)
        -- true anomaly measured from line of nodes -pi..pi
        true
            | mag nodes < epsilon = -- equatorial orbit
                -- measure true anomaly from x-axis
                atan2 y x
            | otherwise =
                -sign (nodes `dot` v) * acos (clamp (-1) 1 (unit r `dot` unit nodes))
            where
            (x, y, _) = r
        arg = 0
    conicTrajectory Parabolic =
        (Orbit p e n mean t0 (i, an, arg), (true, mean, parabolic))
        where
        -- periapsis distance
        q = mu * (1 + cosT) / dot v v
        -- semi-latus rectum
        p = 2 * q
        n = sqrt (mu / (2 * q^3))
        parabolic = trueToEccentricAnomaly e true
        mean = eccentricToMeanAnomaly e parabolic
    conicTrajectory Elliptic =
        (Orbit p e n mean t0 (i, an, arg), (true, mean, eccentric))
        where
        n = sqrt (mu / a^3)
        eccentric = sign true * acos cosE
        mean = eccentricToMeanAnomaly e eccentric
    conicTrajectory Hyperbolic =
        (Orbit p e n mean t0 (i, an, arg), (true, mean, hyperbolic))
        where
        n = sqrt (mu / (-a)^3)
        hyperbolic = sign true * acosh cosE
        mean = eccentricToMeanAnomaly e hyperbolic
    -- semi-major axis
    -- a = 1 / (2 / mag r - dot v v / mu)
    a = p / (1 - e^2)
    -- semi-latus rectum of ellipse or hyperbola
    -- p = a * (1 - e^2)
    p = dot k k / mu

    -- cosine of true anomaly (not well defined for circular orbits!)
    cosT = clamp (-1) 1 $ dot (unit evec) (unit r)
    -- true anomaly -pi..pi
    true = -sign (evec `dot` v) * acos cosT

    -- (hyperbolic) cosine of eccentric anomaly
    cosE = (e + cosT) / (1 + e * cosT)

    -- inclination 0..pi
    i = let (_, _, kz) = k in acos (clamp (-1) 1 (kz / mag k))

    -- longitude of ascending node -pi..pi
    an
        | mag nodes < epsilon = 0     -- equatorial orbit
        | otherwise = atan2 ny nx
        where
        (nx, ny, _) = nodes

    -- argument of periapsis -pi..pi
    arg
        | mag nodes < epsilon =         -- equatorial orbit
            -- measure argument of periapsis from x axis
            atan2 ey ex
        | otherwise =
            sign ez * acos (clamp (-1) 1 (unit nodes `dot` unit evec))
        where
        (ex, ey, ez) = evec
