module Main where

import Test.Framework as TF (defaultMain, testGroup, Test)
import Test.Framework.Providers.QuickCheck2 (testProperty)

import Tests.OrbitTests

tests = [
    testGroup "Orbit tests" [
        testProperty "Not a linear trajectory" notLinear,
        testProperty "Not NaN at time" notNaN,
        testProperty "Not NaN state vectors" notNaN2,
        testProperty "Anomalies within range" anomalyRanges,
        testProperty "Orbital angle element ranges" angleRanges,
        testProperty "Identity" identity,
        testProperty "Gravitational parameter" gravityParameter,
        testProperty "Vis-viva state vectors" visVivaVec,
        testProperty "Vis-viva at time" visVivaTime,
        testProperty "Orbital periodicity" periodicity,
        testProperty "Apsides distance and speed" apsides,
        testProperty "Nodes at equator" nodes,
        testProperty "Velocity test" velocityTest
    ]]

main = defaultMain tests
