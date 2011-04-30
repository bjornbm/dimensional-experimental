module Numeric.Units.Dimensional.Formulae where

import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional.Constants
import qualified Prelude


-- | @dewPoint temperature relativeHumidity@ calculates the dew point
-- based on the August-Roche-Magnus approximation vapor pressure
-- of water in air as a function of temperature. It is considered
-- to have an 1-sigma uncertainty less of 0.4 °C for:
--   0 °C < temperature      < 60 °C,
--   0.01 < relativeHumidity < 1,
--   0 °C < dewPoint t rh    < 50 °C.
--
-- Reference: @http://www.paroscientific.com/dewpoint.htm@
dewPoint :: Floating a => ThermodynamicTemperature a -> Dimensionless a -> ThermodynamicTemperature a
dewPoint t rh = b * gamma t rh / (a - gamma t rh)
  where
    gamma t rh = a * t / (b + t) + log rh
    a = 17.27 *~ one
    b = 237.7 *~ degreeCelsius


-- Arrhenius equation
-- ------------------
-- | @arrheniusEquation e t@ is a variation of Arrhenius equation
-- giving the ratio of collisions exceeding the activation energy @e@
-- (joules per mole) to the total number of collisions in a medium at
-- a temperature @t@.
arrheniusEquation :: Floating a => MolarEnergy a -> ThermodynamicTemperature a -> Dimensionless a
arrheniusEquation e t = exp (negate e / (molarGasConstant * t))

-- | @arrheniusEquation' e t@ is a variation of Arrhenius equation
-- giving the ratio of collisions exceeding the activation energy @e@
-- (joules per molecule) to the total number of collisions in a medium
-- at a temperature @t@.
arrheniusEquation' :: Floating a => Energy a -> ThermodynamicTemperature a -> Dimensionless a
arrheniusEquation' e t = exp (negate e / (boltzmann * t))

-- | Activation energy per mole given two samples of temperature and rate.
activationEnergy :: Floating a
                 => (ThermodynamicTemperature a, Dimensionless a)
                 -> (ThermodynamicTemperature a, Dimensionless a)
                 -> MolarEnergy a
activationEnergy (t1, k1) (t2, k2) = negate molarGasConstant * (log k2 - log k1) / (t2^neg1 - t1^neg1)

-- | Activation energy per molecule given two samples of temperature and rate.
activationEnergy' :: Floating a
                  => (ThermodynamicTemperature a, Dimensionless a)
                  -> (ThermodynamicTemperature a, Dimensionless a)
                  -> Energy a
activationEnergy' (t1, k1) (t2, k2) = negate boltzmann * (log k2 - log k1) / (t2^neg1 - t1^neg1)

-- | @http://en.wikipedia.org/wiki/Calorie@
thermochemicalCalorie :: Fractional a => Unit DEnergy a
thermochemicalCalorie = prefix 4.184 joule
