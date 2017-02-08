
-- The convention in this module is that a @C@ denotes cartesian coordinates and an @S@ denotes spherical coordinates.

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Numeric.Units.Dimensional.LinearAlgebra.PosVel where

import qualified Prelude
import Numeric.NumType.DK.Integers (TypeInt (Pos2))
import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional.LinearAlgebra
import Numeric.Units.Dimensional.LinearAlgebra.VectorAD (applyLinear)


-- Type synonyms for clearer documentation.

-- | Angle from Z-axis.
type DZenith = DPlaneAngle; type Zenith = Angle
-- | Angle from X-axis towards Y-axis (positive about Z).
type DRightAscension = DPlaneAngle; type RightAscension = Angle


-- Some type synonyms for convenience.

type Vec3 d = Vec d 3
type Car d = Vec3 d  -- ^ Cartesian vector (x,y,z).
type CPos = Car DLength
type CVel = Car DVelocity

-- TODO is this general spherical coordinate useful or should I jusst specialize to Radius?
-- ^ Spherical coordinates.
data Sph d a = Sph { magnitude      :: Quantity d a
                   , zenith         :: PlaneAngle a
                   , rightAscension :: PlaneAngle a
                   } deriving (Eq)
deriving instance (KnownDimension d, Show a, Fractional a) => Show (Sph d a)  -- Needed since d unknown.
type SPos = Sph DRadius
data SVel a = SVel (Velocity a) (AngularVelocity a) (AngularVelocity a)
  deriving (Eq, Show)



-- Querying
-- --------
--
-- Cartesian coordinates.

x :: Car d a -> Quantity d a
x = vElemAt n0
y :: Car d a -> Quantity d a
y = vElemAt n1
z :: Car d a -> Quantity d a
z = vElemAt n2

-- Spherical coordinates.
colatitude  :: Sph d a -> Zenith a
colatitude  = zenith
polarAngle  :: Sph d a -> Zenith a
polarAngle  = zenith
latitude    :: Floating a => Sph d a -> Angle a
latitude s  = tau / _4 - colatitude s
declination :: Floating a => Sph d a -> Angle a
declination = latitude
longitude :: Sph d a -> RightAscension a
longitude = rightAscension
hourAngle :: Sph d a -> RightAscension a
hourAngle = rightAscension


-- Converting
-- ----------
-- Converts a cartesian position vector into a spherical position vector.

c2s :: RealFloat a => Car d a -> Sph d a
c2s c = Sph r zen az
  where
    (x, y, z) = toTuple c
    r   = norm c  -- sqrt (x^pos2 + y^pos2 + z^pos2)
    zen = if r == _0 then _0 else acos (z / r)
    az  = atan2 y x


-- Converts a spherical position vector into a cartesian position vector.

s2c :: Floating a => Sph d a -> Car d a
s2c s = x <: y <:. z
  where
    Sph r zen az = s
    x = r * sin zen * cos az
    y = r * sin zen * sin az
    z = r * cos zen

-- Convert declination/latitude to zenith.

toZenith :: Floating a => Angle a -> Zenith a
toZenith decl = tau / _4 - decl


-- Position vectors
-- ================

type DRadius = DLength; type Radius = Length

radius :: SPos a -> Radius a
radius = magnitude


-- Position and velocity ephemeris
-- ===============================
-- Data type combining position and velocity into a state vector (minus epoch).

type CPosVel a = (CPos a, CVel a)
type SPosVel a = (SPos a, SVel a)


-- TODO the below doesn't work due to SPos and SVel not being representable as homogeneous vectors. :(
-- Conversions
-- -----------

c2sEphem :: RealFloat a => CPosVel a -> SPosVel a
c2sEphem = pvf . applyLinear (fp . c2s) -- unlinearize (c2s . linearize c :: RealFloat b => Time b -> SPos b)
  where
    fp (Sph r zen az) = r / (1 *~ meter) <: zen <:. az
    pvf (p, v) = (pf p, vf v)
      where
        pf v = Sph  (r  * (1 *~ meter)) zen  az  where [r , zen , az ] = listElems v
        vf v = SVel (r' * (1 *~ meter)) zen' az' where [r', zen', az'] = listElems v

s2cEphem :: Floating a => SPosVel a -> CPosVel a
s2cEphem = applyLinear (s2c . pf) . fpv -- unlinearize (s2c . linearize s :: RealFloat b => Time b -> CPos b)
  where
    pf v = Sph (r * (1 *~ meter)) zen az where [r, zen, az] = listElems v
    fpv (p,v) = (fp p, fv v)
      where
        fp (Sph r zen az) = r / (1 *~ meter) <: zen <:. az
        fv (SVel r' zen' az') = r' / (1 *~ meter) <: zen' <:. az'
