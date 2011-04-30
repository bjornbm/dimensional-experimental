{-# LANGUAGE TypeOperators #-}

module Numeric.Units.Dimensional.LinearAlgebra.Rotation where

import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional.LinearAlgebra
import qualified Prelude

-- | The type of homogenous vectors with three elements.
type Homo3 d = Vec (d:*:d:*.d)

-- | Cartesian unit vectors.
unit_x, unit_y, unit_z :: Num a => Homo3 DOne a
unit_x = vCons _1 $ vCons _0 $ vSing _0
unit_y = vCons _0 $ vCons _1 $ vSing _0
unit_z = vCons _0 $ vCons _0 $ vSing _1


-- Rotation matrices (cartesian)
-- =============================

-- | Convenience type for homogeneous 3x3 matrices.
type Homo33 d = Mat ((d:*:d:*.d) :*:
                     (d:*:d:*.d) :*.
                     (d:*:d:*.d))

-- Rotation matrices. Rotates a vector by the given angle (analogous
-- to rotating the coordinate system in opposite direction).

rotX :: Floating a => PlaneAngle a -> Homo33 DOne a
rotX a = consRow   (vCons _1 $ vCons _0      $ vSing _0)
       $ consRow   (vCons _0 $ vCons (cos a) $ vSing (negate (sin a)))
       $ rowMatrix (vCons _0 $ vCons (sin a) $ vSing (cos a))

rotY :: Floating a => PlaneAngle a -> Homo33 DOne a
rotY a = consRow   (vCons (cos a)          $ vCons _0 $ vSing (sin a))
       $ consRow   (vCons _0               $ vCons _1 $ vSing _0)
       $ rowMatrix (vCons (negate (sin a)) $ vCons _0 $ vSing (cos a))

rotZ :: Floating a => PlaneAngle a -> Homo33 DOne a
rotZ a = consRow   (vCons (cos a) $ vCons (negate (sin a)) $ vSing _0)
       $ consRow   (vCons (sin a) $ vCons (cos a)          $ vSing _0)
       $ rowMatrix (vCons _0      $ vCons _0               $ vSing _1)
