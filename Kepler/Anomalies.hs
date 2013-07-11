module Kepler.Anomalies (
    meanToEccentricAnomaly,
    meanToTrueAnomaly,
    eccentricToTrueAnomaly,
    eccentricToMeanAnomaly,
    trueToEccentricAnomaly,
    trueToMeanAnomaly,
    meanAnomalyAtTime,
    anomalies
    ) where

import Kepler.MathOps (epsilon, sign, angleClamp)
import Kepler.Orbit (Orbit(..), OrbitType(..), meanAnomalyAtEpoch, meanMotion, eccentricType, orbitClosed)

meanToEccentricAnomaly e mean =
    anomaly (eccentricType e)
    where
    anomaly Circular =
        mean
    anomaly Parabolic =
        -- D^3 / 3 + D - M = 0
        -- the only (always) real-valued root of D is: [wolfram alpha]
        (num / denom)**(1/3) - (denom / num)**(1/3)
        where
        num = 2
        denom = sqrt (9 * mean^2 + 4) - 3 * mean
    anomaly _ =
        iterate f mean steps
        where
        steps
            | e < 0.3 = 10
            | e < 0.9 = 10
            | e < 1 = 20
            | e > 1 = 30
        f x
            | e < 0.3 =
                mean + e * sin x
            | e < 0.9 =
                x + (mean + e * sin x - x) / (1 - e * cos x )
            | e < 1 =
                let s = e * sin x
                    c = e * cos x
                    f0 = x - s - mean
                    f1 = 1 - c
                    f2 = s
                in
                    x + (-5) * f0 / (f1 + sign(f1) * sqrt(abs(16 * f1 * f1 - 20 * f0 * f2)))
            | e > 1 =
                let s = e * sinh x
                    c = e * cosh x
                    f0 = s - x - mean
                    f1 = c - 1
                    f2 = s
                in
                    x + (-5) * f0 / (f1 + sign(f1) * sqrt(abs(16 * f1 * f1 - 20 * f0 * f2)))
        sign x
            | x < 0 = -1
            | otherwise = 1
        iterate _ x 0 = x
        iterate f x n =
            if abs (eccentricToMeanAnomaly e x - mean) < epsilon then x else iterate f (f x) (n-1)

eccentricToMeanAnomaly e eccentric =
    case (eccentricType e) of
        Circular ->
            eccentric
        Elliptic ->
            eccentric - e * sin eccentric
        Parabolic ->
            eccentric^3 / 3 + eccentric
        Hyperbolic ->
            e * sinh eccentric - eccentric

eccentricToTrueAnomaly e eccentric =
    case (eccentricType e) of
        Circular ->
            eccentric
        Elliptic ->
            atan2 (sqrt (1 - e^2) * sin eccentric) (cos eccentric - e)
        Parabolic ->
            2 * atan eccentric
        Hyperbolic ->
            2 * atan (sqrt ((e + 1) / (e - 1)) * tanh (eccentric/2) )

trueToEccentricAnomaly e true =
    case (eccentricType e) of
        Circular ->
            true
        Elliptic ->
            atan2 (sqrt (1 - e^2) * sin true) (cos true + e)
        Parabolic ->
            atan (true / 2)
        Hyperbolic ->
            sign true * acosh ((e + cos true) / (1 + e * cos true))

meanToTrueAnomaly e = eccentricToTrueAnomaly e . meanToEccentricAnomaly e
trueToMeanAnomaly e = eccentricToMeanAnomaly e . trueToEccentricAnomaly e

meanAnomalyAtTime orbit t =
    if orbitClosed orbit then angleClamp m else m
    where
    m = meanAnomalyAtEpoch orbit + (t - orbitEpoch orbit) * meanMotion orbit

anomalies elems@(Orbit _ e _ _ _ _) t =
    (mean, eccentric, true)
    where
    mean = meanAnomalyAtTime elems t
    eccentric = meanToEccentricAnomaly e mean
    true = eccentricToTrueAnomaly e eccentric
