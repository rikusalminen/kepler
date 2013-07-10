module Tests.InterceptTests where

import Test.QuickCheck
import Debug.Trace

import Kepler.MathOps
import Kepler.Orbit
import Kepler.OrbitState
import Kepler.OrbitIntercept

import Tests.Arbitrary

interceptTest (TestIntercept (o1, o2) (t0, t1) t) =
    err < epsilon
    where
    intercepts = []
    ((r1, _), (r2, _)) = (orbitStateAtTime o1 t, orbitStateAtTime o2 t)
    dr = sub r2 r1
    err = (dot dr dr) / (mag r1 + mag r2)^2
    --intercepts = orbitIntercept (o1, o2) (t0, t1) 0
