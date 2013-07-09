module Main where

import Kepler.Orbit
import Kepler.Anomalies
import Kepler.OrbitState

import Control.Monad (forM_)

main = do
    forM_ states write
    where
    mu = 40.0e13
    testElems = [
        --orbitStateToElements mu 0.0 (q, 0.0, 0.0) (0.0, 2.0 * vEsc, 0.0) 0.0,
        orbitStateToElements mu ((q, 0.0, 0.0), (0.0, 1.1 * vEsc, 0.0)) 0.0,
        orbitStateToElements mu ((q, 0.0, 0.0), (0.0, vEsc, 0.0)) 0.0,
        orbitStateToElements mu ((q, 0.0, 0.0), (0.0, 1.1 * vCirc, 0.0)) 0.0,
        orbitStateToElements mu ((q, 0.0, 0.0), (0.0, vCirc, 0.0)) 0.0,
        orbitStateToElements mu ((q, 0.0, 0.0), (0.0, 0.0, vCirc)) 0.0,
        orbitStateToElements mu ((0.0, q, 0.0), (0.0, 0.0, vCirc)) 0.0,
        orbitStateToElements mu ((q, 0.0, 0.0), (0.0, 0.75 * vCirc, 0.0)) 0.0
        ]
    q = 7500.0e3
    vEsc = sqrt (2.0 * mu / q)
    vCirc = sqrt (mu / q)
    num = 64
    write ((x, y, z), (vx, vy, vz)) = do
        putStrLn $ show x ++ "\t" ++ show y ++ "\t" ++ show z ++ "\t" ++ show vx ++ "\t" ++ show vy ++ "\t" ++ show vz
    states = concat [map (state elems) times | (elems, _) <- testElems]
    times = [(m / num) | m <- [0..num]]
    state orbit t' =
        orbitState orbit anoms
        where
        t = (if orbitClosed orbit then t'*period else (t'*period - period / 2))
        anoms = anomalies orbit (orbitEpoch orbit + t)
        period = if orbitClosed orbit then orbitalPeriod orbit else 60*60 * 24
