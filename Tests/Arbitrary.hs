module Tests.Arbitrary where

import Test.QuickCheck

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
            (40, choose (1, 5))
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
            (40, choose (1.0, 5.0))
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
        let t = orbitEpoch orbit + m / meanMotion orbit

        n <- choose (-1000, 1000 :: Int)
        let time
                | orbitClosed orbit =
                    t + (fromIntegral n) * orbitalPeriod orbit
                | otherwise =
                    t

        return $ TestTime orbit time
