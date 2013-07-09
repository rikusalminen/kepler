module Kepler.Orbit (
    Orbit(..),
    OrbitType(..),
    orbitType,
    eccentricType,
    semiMajorAxis, semiMinorAxis,
    periapsis, apoapsis,
    orbitGravitationalParameter,
    periapsisVel, apoapsisVel,
    orbitClosed,
    orbitalPeriod,
    orbitMatrix,
    orientationMatrix
    ) where

import Kepler.MathOps (epsilon, inf)

data Orbit = Orbit {
    semiLatusRectum :: Double,
    eccentricity :: Double,
    meanMotion :: Double,
    meanAnomalyAtEpoch :: Double,
    orbitEpoch :: Double,
    orbitOrientation :: (Double, Double, Double)
    }
    deriving Show

data OrbitType = Circular | Elliptic | Parabolic | Hyperbolic
    deriving (Show, Eq, Ord)

orbitType = eccentricType . eccentricity

eccentricType e
    | e < epsilon = Circular
    | abs (e - 1) < epsilon = Parabolic
    | e < 1 = Elliptic
    | e > 1 = Hyperbolic

semiMajorAxis (Orbit p e _ _ _ _) =
    case eccentricType e of
        Parabolic ->
            inf
        _ ->
            p / (1 - e^2)

semiMinorAxis (Orbit p e _ _ _ _) =
    case eccentricType e of
        Parabolic ->
            inf
        Hyperbolic ->
            -p / sqrt (e^2 - 1)
        _ ->
            p / sqrt (1 - e^2)

periapsis (Orbit p e _ _ _ _) =
    p / (1 + e)

apoapsis orbit@(Orbit p e _ _ _ _)
    | orbitClosed orbit =
        p / (1 - e)
    | otherwise = inf

orbitGravitationalParameter (Orbit p e n _ _ _) =
    case eccentricType e of
        Parabolic ->
            n^2 * p^3 / 4
        _ ->
            n^2 * (p / abs (1 - e^2))^3

periapsisVel orbit@(Orbit p e _ _ _ _) =
    sqrt ((mu / p) * (1 + e)^2)
    where
    mu = orbitGravitationalParameter orbit

apoapsisVel orbit@(Orbit p e _ _ _ _)
    | orbitClosed orbit =
        sqrt (mu / p * (1 - e)^2)
    | otherwise =
        inf
    where
    mu = orbitGravitationalParameter orbit

orbitClosed orbit = 1 - eccentricity orbit > epsilon

orbitalPeriod (Orbit _ e n _ _ _)
    | e < 1 = 2 * pi / n
    | otherwise = inf

orientationMatrix (i, an, arg) =
    (
        (
            (cos arg * cos an) - (sin arg * sin an * cos i),
            -(cos arg * sin an * cos i) - (sin arg * cos an),
            (sin an * sin i)
        ),
        (
            (sin arg * cos an * cos i) + (cos arg * sin an),
            (cos arg * cos an * cos i) - (sin arg * sin an),
            -(cos an * sin i)
        ),
        (
            sin arg * sin i,
            cos arg * sin i,
            cos i
        )
    )

orbitMatrix = orientationMatrix . orbitOrientation
