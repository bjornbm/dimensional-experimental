
The convention in this module is that a @C@ denotes cartesian coordinates and an @S@ denotes spherical coordinates.

> {-# OPTIONS_GHC -fglasgow-exts #-}
> {-# LANGUAGE Haskell98 #-}

> module Numeric.Units.Dimensional.LinearAlgebra.PosVel where

> import qualified Prelude
> import Numeric.NumType (Pos2)
> import Numeric.Units.Dimensional.Prelude
> import Numeric.Units.Dimensional.LinearAlgebra
> import Numeric.Units.Dimensional.LinearAlgebra.VectorAD (applyLinear)


Type synonyms for clearer documentation.

> -- | Angle from Z-axis.
> type DZenith = DPlaneAngle; type Zenith = Angle
> -- | Angle from X-axis towards Y-axis (positive about Z).
> type DRightAscension = DPlaneAngle; type RightAscension = Angle


Some type synonyms for convenience.

> type Vec3 d1 d2 d3 = Vec (d1:*:d2:*.d3)
> type Car d = Vec3 d d d  -- ^ Cartesian vector (x,y,z).
> type Sph d = Vec3 d DPlaneAngle DPlaneAngle  -- ^ Spherical vector.
> type CPos = Car DLength
> type CVel = Car DVelocity
> type SPos = Sph DRadius
> type SVel = Vec3 DVelocity DAngularVelocity DAngularVelocity



Querying
--------

Cartesian coordinates.

> x :: Car d a -> Quantity d a
> x = vElemAt zero
> y :: Car d a -> Quantity d a
> y = vElemAt pos1
> z :: Car d a -> Quantity d a
> z = vElemAt pos2

Spherical coordinates.

> magnitude :: Sph d a -> Quantity d a
> magnitude = vElemAt zero

> zenith :: Sph d a -> Zenith a
> zenith = vElemAt pos1
> colatitude  :: Sph d a -> Zenith a
> colatitude  = zenith
> polarAngle  :: Sph d a -> Zenith a
> polarAngle  = zenith
> latitude    :: Floating a => Sph d a -> Angle a
> latitude s  = pi / _2 - colatitude s
> declination :: Floating a => Sph d a -> Angle a
> declination = latitude

> rightAscension :: Sph d a -> RightAscension a
> rightAscension = vElemAt pos2
> longitude      :: Sph d a -> RightAscension a
> longitude      = rightAscension
> hourAngle      :: Sph d a -> RightAscension a
> hourAngle      = rightAscension


Converting
----------
Converts a cartesian position vector into a spherical position vector.

> c2s :: (Div d d DOne, Mul d d d2, Root d2 Pos2 d, RealFloat a)
>     => Car d a -> Sph d a
> c2s c = fromTuple (r, zen, az)
>   where
>     (x, y, z) = toTuple c
>     r   = vNorm c  -- sqrt (x^pos2 + y^pos2 + z^pos2)
>     zen = if r == _0 then _0 else acos (z / r)
>     az  = atan2 y x


Converts a spherical position vector into a cartesian position vector.

> s2c :: (Mul d DOne d) => Floating a => Sph d a -> Car d a
> s2c s = fromTuple (x, y, z)
>   where
>     (r, zen, az) = toTuple s
>     x = r * sin zen * cos az
>     y = r * sin zen * sin az
>     z = r * cos zen


Position vectors
================

> type DRadius = DLength; type Radius = Length

> radius :: SPos a -> Radius a
> radius = magnitude


Position and velocity ephemeris
===============================
Data type combining position and velocity into a state vector (minus epoch).

> type CPosVel a = (CPos a, CVel a)
> type SPosVel a = (SPos a, SVel a)


Conversions
-----------

> c2sEphem :: RealFloat a => CPosVel a -> SPosVel a
> c2sEphem = applyLinear c2s -- unlinearize (c2s . linearize c :: RealFloat b => Time b -> SPos b)

> s2cEphem :: RealFloat a => SPosVel a -> CPosVel a
> s2cEphem = applyLinear s2c -- unlinearize (s2c . linearize s :: RealFloat b => Time b -> CPos b)


