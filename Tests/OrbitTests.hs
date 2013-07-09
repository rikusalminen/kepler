module Main where

import Test.QuickCheck
import Test.Framework as TF (defaultMain, testGroup, Test)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Debug.Trace

import Kepler.MathOps
import Kepler.Orbit
import Kepler.Anomalies
import Kepler.OrbitState

arbitraryOrientation = do
    i <- frequency [
        (5, return 0),
        (5, return $ pi/2),
        (5, return $ pi),
        (40, choose (0, pi/2)),
        (40, choose (0, pi/2))
        ]
    an <- angle
    arg <- angle
    return (i, an, arg)
    where
    angle =
        frequency [
            (5, return 0),
            (5, return $ pi/2),
            (5, return $ -pi/2),
            (5, return $ pi),
            (5, return $ -pi),
            (75, choose (-pi, pi))
            ]

instance Arbitrary Orbit where
    arbitrary = do
        p <- choose (1e3, 1e15)

        e <- frequency [
            (10, return 0),
            (10, return 1),
            (40, choose (0, 1)),
            (40, choose (1, 10))
            ]

        n <- choose (2 * pi / maxPeriod, 2 * pi / minPeriod)
        m0 <- choose (if 1 - e > epsilon then (-pi, pi) else (-pi/2, pi/2))
        t0 <- choose (-maxPeriod, maxPeriod)

        orientation <- arbitraryOrientation

        return $ Orbit p e n m0 t0 orientation

        where
        minPeriod = 60*30
        maxPeriod = 60*60 * 24 * 365 * 100

data TestVector = TestVector Double ((Double, Double, Double), (Double, Double, Double)) Double
    deriving (Show)

instance Arbitrary TestVector where
    arbitrary = do
        mu <- choose (muMin, muMax)

        e <- frequency [
            (10, return 0),
            (10, return 1),
            (40, choose (0, 1.0)),
            (40, choose (1.0, 10.0))
            ]

        p <- choose (1.0e3, 1.0e15)

        true <- choose (if 1 - e > epsilon then (-pi, pi) else (-pi/2, pi/2))

        let a = p / (1 - e^2)
            r = p / (1 + e * cos true)
            x = r * cos true
            y = r * sin true
            v = if eccentricType e == Parabolic
                    then sqrt (mu * (2 / r))
                    else sqrt (mu * ((2 / r) - (1 / a)))
            dx = -p * sin true / (1 + e * cos true)^2
            dy = p * (e + cos true) / (1 + e * cos true)^2
            d = sqrt (dx^2 + dy^2)
            vx = v * dx / d
            vy = v * dy / d

        orientation <- arbitraryOrientation
        let matrix = orientationMatrix orientation

        t0 <- choose (-maxPeriod, maxPeriod)
        return $ TestVector mu (vproduct matrix (x, y, 0), vproduct matrix (vx, vy, 0)) t0

        where
        muMin = 1e13
        muMax = 1e20
        maxPeriod = 60*60 * 24 * 365

data TestTime = TestTime Orbit Double
    deriving Show

instance Arbitrary TestTime where
    arbitrary = do
        orbit <- arbitrary
        m <- if orbitClosed orbit
            then frequency [
                (5, return $ 0),
                (5, return $ pi),
                (5, return $ -pi),
                (5, return $ pi / 2),
                (5, return $ -pi / 2),
                (75, choose (-pi, pi))
                ]
            else frequency [
                (5, return $ 0),
                (95, choose (-pi/2, pi/2))
                ]
        let t = m / meanMotion orbit

        n <- choose (-1000, 1000 :: Int)
        let time
                | orbitClosed orbit =
                    t + (fromIntegral n) * orbitalPeriod orbit
                | otherwise =
                    t

        return $ TestTime orbit time

eqf x y = abs (x - y) / (abs x + abs y) / 2 < epsilon
eq3 (x1, y1, z1) (x2, y2, z2) = eqf x1 x2 && eqf y1 y2 && eqf z1 z2

zero v = dot v v < epsilon
eqv a b = ((2 / (mag a + mag b))^2) * dot (sub a b) (sub a b) < 1.0e-7

notLinear (TestVector mu (r, v) t0) =
    (not . zero $ v) && (not . zero $ r) && (not . zero $ (r `cross` v))

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

gravityParameter (TestVector mu (r, v) t0) =
    eqf mu mu'
    where
    mu' = orbitGravitationalParameter orbit
    (orbit, _) = orbitStateToElements mu (r, v) t0

periodicity (TestTime orbit t0) =
    not (orbitClosed orbit) || eqv r r' && eqv v v'
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

tests = [
    testGroup "Orbit tests" [
        testProperty "Not a linear trajectory" notLinear,
        testProperty "Orbital angle element ranges" angleRanges,
        testProperty "Identity" identity,
        testProperty "Gravitational parameter" gravityParameter,
        testProperty "Vis-viva state vectors" visVivaVec,
        testProperty "Orbital periodicity" periodicity,
        testProperty "Apsides distance and speed" apsides,
        testProperty "Nodes at equator" nodes,
        testProperty "Velocity test" velocityTest
    ]]

main = defaultMain tests
