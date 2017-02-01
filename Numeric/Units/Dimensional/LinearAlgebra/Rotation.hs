{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE DataKinds #-}

module Numeric.Units.Dimensional.LinearAlgebra.Rotation where

import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional.LinearAlgebra
import qualified Prelude

-- | The type of homogenous vectors with three elements.
type Homo3 d = Vec d 3

-- | Cartesian unit vectors.
unit_x, unit_y, unit_z :: Num a => Homo3 DOne a
unit_x = vCons _1 $ vCons _0 $ vSing _0
unit_y = vCons _0 $ vCons _1 $ vSing _0
unit_z = vCons _0 $ vCons _0 $ vSing _1


-- Rotation matrices (cartesian)
-- =============================

-- | Convenience type for homogeneous 3x3 matrices.
type Homo33 d = Mat d 3 3

-- Rotation matrices. Rotates a vector by the given angle (analogous
-- to rotating the coordinate system in opposite direction).

rotX :: Floating a => PlaneAngle a -> Homo33 DOne a
rotX a =   (_1 <:    _0 <:.            _0 )
       |:  (_0 <: cos a <:. negate (sin a))
       |:. (_0 <: sin a <:.         cos a )

rotY :: Floating a => PlaneAngle a -> Homo33 DOne a
rotY a =   (        cos a  <: _0 <:. sin a)
       |:  (           _0  <: _1 <:.    _0)
       |:. (negate (sin a) <: _0 <:. cos a)

rotZ :: Floating a => PlaneAngle a -> Homo33 DOne a
rotZ a =   (cos a <: negate (sin a) <:. _0)
       |:  (sin a <:         cos a  <:. _0)
       |:. (   _0 <:            _0  <:. _1)

-- | @rotV v a@ generates the rotation matrix for a counterclockwise
-- rotation through the angle @a@ about an axis in R^3, which is
-- determined by the arbitrary unit vector @v@ that has its initial
-- point at the origin (from Anton, Rorres, 1994, pages 190-191).
rotV :: Floating a => Homo3 DOne a -> Angle a -> Homo33 DOne a
rotV v t =   (a*a*omct +   cos t <: a*b*omct - c*sin t <:. a*c*omct + b*sin t)
         |:  (b*a*omct + c*sin t <: b*b*omct +   cos t <:. b*c*omct - a*sin t)
         |:. (c*a*omct - b*sin t <: c*b*omct + a*sin t <:. c*c*omct +   cos t)
  where
    (a,b,c) = toTuple v
    omct    = _1 - cos t
