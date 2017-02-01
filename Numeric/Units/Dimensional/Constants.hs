{-# LANGUAGE DataKinds #-}

module Numeric.Units.Dimensional.Constants where

import Numeric.Units.Dimensional.Prelude
import Numeric.NumType.DK.Integers (TypeInt (Neg4, Neg3, Neg2, Neg1, Zero, Pos1, Pos2, Pos3))
import qualified Prelude

-- Dim l m t i th n j

-- CODATA 2010

-- Universal constants
-- ===================

-- | Magnetic constant (mu_0).
permeability :: Floating a => Permeability a
permeability = 4e-7 *~ (volt * second / ampere / meter) * pi  -- in vacuum.
magneticConstant :: Floating a => Permeability a
magneticConstant = permeability  -- CODATA name.

-- | Electric constant (epsilon_0).
permittivity :: Floating a => Permittivity a
permittivity = permeability^neg1 * speedOfLight^neg2  -- in vacuum.
electricConstant :: Floating a => Permittivity a
electricConstant = permittivity  -- CODATA name.

-- | Newtonian constant of gravitation (G).
gravitationalConstant :: Fractional a => Quantity (Dim Pos3 Neg1 Neg2 Zero Zero Zero Zero) a
gravitationalConstant = 6.67384e-11 *~ (newton * meter^pos2 / kilo gram^pos2)
-- Uncertainty:              80

-- | Planck constant (h).
planck :: Fractional a => Quantity (Dim Pos2 Pos1 Neg1 Zero Zero Zero Zero) a
planck = 6.62606957e-34 *~ (joule * second)
-- Uncertainty:  29

-- | Planck length (l_P).
planckLength :: Fractional a => Length a
planckLength = 1.616199e-35 *~ meter

-- | Planck mass (m_P).
planckMass :: Fractional a => Mass a
planckMass = 2.17651e-8  *~ kilo gram
-- Uncertainty:   13

-- | Planck temperature (T_P).
planckTemperature :: Fractional a => ThermodynamicTemperature a
planckTemperature = 1.416833e-32 *~ kelvin
-- Uncertainty:           85

-- | Planck time (t_P).
planckTime :: Fractional a => Time a
planckTime = 5.39106e-44 *~ second
-- Uncertainty:   32

-- | Spead of light in vacuum (c, c_0).
speedOfLight :: Fractional a => Velocity a
speedOfLight = 299792458 *~ (meter / second) -- in vacuum.


-- Atomic and nuclear constants
-- ============================

-- | Electron mass (m_e).
electronMass :: Fractional a => Mass a
electronMass = 9.10938291e-31 *~ kilo gram
-- Uncertainty:        40

-- | Elementary charge (e).
elementaryCharge :: Fractional a => ElectricCharge a
elementaryCharge = 1.602176565e-19 *~ coulomb
-- Uncertainty:             35

-- | Rydberg constant (R_inf).
-- TODO dimensional-codata?
rydberg :: Floating a => WaveNumber a
rydberg = electronMass * elementaryCharge^pos4 / (_8 * permittivity^pos2 * planck^pos3 * speedOfLight)
-- Using Double: 1.0973731593928682e7 m^-1
-- Per CODATA:   1.0973731568539(55)  m^-1

-- Physico-chemical constants
-- ==========================

-- Dim l m t i th n j

-- | Avogadro constant (N_A, L).
avogadro :: Fractional a => Quantity (Dim Zero Zero Zero Zero Zero Neg1 Zero) a
avogadro = 6.02214129e23 *~ mole^neg1
-- Uncertainty:    27

-- | Boltzmann constant (k).
boltzmann :: Fractional a => Entropy a
boltzmann = 1.3806488e-23 *~ (joule / kelvin)
-- Uncertainty:    13

-- | Molar gas constant (R = k * N_A).
-- TODO dimensional-codata?
molarGasConstant :: Fractional a => MolarEntropy a
molarGasConstant = boltzmann * avogadro

-- | Stefan-Boltzmann constant (sigma).
-- TODO dimensional-codata?
stefanBoltzmann :: Floating a => Quantity (Dim Zero Pos1 Neg3 Zero Neg4 Zero Zero) a
stefanBoltzmann = _2 * pi ^ pos5 * boltzmann^pos4 / (15*~one) / planck^pos3 / speedOfLight ^pos2
