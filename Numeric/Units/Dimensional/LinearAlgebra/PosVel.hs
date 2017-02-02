
-- The convention in this module is that a @C@ denotes cartesian coordinates and an @S@ denotes spherical coordinates.

{-# OPTIONS_GHC -fglasgow-exts #-}
{-# LANGUAGE Haskell98 #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}

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
data Sph d a = Sph (Quantity d a) (PlaneAngle a) (PlaneAngle a)  -- ^ Spherical coordinates.
type CPos = Car DLength
type CVel = Car DVelocity
type SPos = Sph DRadius
data SVel a = SVel (Velocity a) (AngularVelocity a) (AngularVelocity a)



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

magnitude :: Sph d a -> Quantity d a
magnitude (Sph r _ _)= r

zenith :: Sph d a -> Zenith a
zenith (Sph _ z _)= z
colatitude  :: Sph d a -> Zenith a
colatitude  = zenith
polarAngle  :: Sph d a -> Zenith a
polarAngle  = zenith
latitude    :: Floating a => Sph d a -> Angle a
latitude s  = tau / _4 - colatitude s
declination :: Floating a => Sph d a -> Angle a
declination = latitude

rightAscension :: Sph d a -> RightAscension a
rightAscension (Sph _ _ ra)= ra
longitude      :: Sph d a -> RightAscension a
longitude      = rightAscension
hourAngle      :: Sph d a -> RightAscension a
hourAngle      = rightAscension


-- Converting
-- ----------
-- Converts a cartesian position vector into a spherical position vector.

c2s :: RealFloat a => Car d a -> Sph d a
c2s c = Sph r zen az
  where
    (x, y, z) = toTuple c
    r   = vNorm c  -- sqrt (x^pos2 + y^pos2 + z^pos2)
    zen = if r == _0 then _0 else acos (z / r)
    az  = atan2 y x


-- Converts a spherical position vector into a cartesian position vector.

s2c :: Floating a => Sph d a -> Car d a
s2c s = fromTuple (x, y, z)
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

c2sEphem :: forall a . RealFloat a => CPosVel a -> SPosVel a
c2sEphem = pvf . applyLinear (fp . c2s) -- unlinearize (c2s . linearize c :: RealFloat b => Time b -> SPos b)
  where
    fp :: RealFloat c => SPos c -> Vec3 DOne c
    fp (Sph r zen az) = r / (1 *~ meter) <: zen <:. az
    pvf :: (Vec3 DOne a, Vec3 DFrequency a) -> (SPos a, SVel a)
    pvf (p, v) = (pf p, vf v)
      where
        pf :: RealFloat c => Vec3 DOne c -> SPos c
        pf v = Sph (r * (1 *~ meter)) zen az where [r, zen, az] = toList v
        vf :: RealFloat c => Vec3 DFrequency c -> SVel c
        vf v = SVel (r' * (1 *~ meter)) zen' az' where [r', zen', az'] = toList v

s2cEphem :: forall a . RealFloat a => SPosVel a -> CPosVel a
s2cEphem = applyLinear (s2c . pf) . fpv -- unlinearize (s2c . linearize s :: RealFloat b => Time b -> CPos b)
  where
    pf :: RealFloat c => Vec3 DOne c -> SPos c
    pf v = Sph (r * (1 *~ meter)) zen az where [r, zen, az] = toList v
    fpv :: (SPos a, SVel a) -> (Vec3 DOne a, Vec3 DFrequency a)
    fpv (p,v) = (fp p, fv v)
      where
        fp :: SPos a -> Vec3 DOne a
        fp (Sph r zen az) = r / (1 *~ meter) <: zen <:. az
        fv :: SVel a -> Vec3 DFrequency a
        fv (SVel r' zen' az') = r' / (1 *~ meter) <: zen' <:. az'
