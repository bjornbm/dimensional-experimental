module Numeric.Units.Dimensional.Constants where

import Numeric.Units.Dimensional.Prelude
import Numeric.NumType (Neg3, Neg2, Neg1, Zero, Pos1, Pos2, Pos3)
import qualified Prelude

-- Dim l m t i th n j

-- Physics handbook p12
-- Empty Space
speedOfLight :: Fractional a => Velocity a
speedOfLight = 2.99792458e8 *~ (meter / second) -- in vacuum.
permeability :: Floating a => Permeability a
permeability = 4e-7 *~ (volt * second / ampere / meter) * pi
permittivity :: Floating a => Permittivity a
permittivity = permeability^neg1 * speedOfLight^neg2

-- Gravitation
gravitationalConstant :: Fractional a => Quantity (Dim Pos3 Neg1 Neg2 Zero Zero Zero Zero) a
gravitationalConstant = 6.67259e-11 *~ (newton * meter^pos2 / kilo gram^pos2)

-- Particle Masses
electronMass :: Fractional a => Mass a
electronMass = 9.109390e-31 *~ kilo gram

-- Physics handbook p13
elementaryCharge :: Fractional a => ElectricCharge a
elementaryCharge = 1.6021773e-19 *~ coulomb

-- Physics handbook p14
planck :: Fractional a => Quantity (Dim Pos2 Pos1 Neg1 Zero Zero Zero Zero) a
planck = 6.626076e-34 *~ (joule * second)
rydberg :: Floating a => WaveNumber a
rydberg = electronMass * elementaryCharge^pos4 / (_8 * permittivity^pos2 * planck^pos3 * speedOfLight)
boltzmann :: Fractional a => Entropy a
boltzmann = 1.38066e-23 *~ (joule / kelvin)

-- Dim l m t i th n j

-- Physics handbook p15
planckLength :: Fractional a => Length a
planckLength      = 1.616e-35 *~ meter
planckMass :: Fractional a => Mass a
planckMass        = 2.177e-8  *~ kilo gram
planckTime :: Fractional a => Time a
planckTime        = 5.391e-44 *~ second
planckTemperature :: Fractional a => ThermodynamicTemperature a
planckTemperature = 1.417e-32 *~ kelvin

avogrado :: Fractional a => Quantity (Dim Zero Zero Zero Zero Zero Neg1 Zero) a
avogrado = 6.022137e23 *~ mole^neg1
molarGasConstant :: Fractional a => MolarEntropy a
molarGasConstant = boltzmann * avogrado
