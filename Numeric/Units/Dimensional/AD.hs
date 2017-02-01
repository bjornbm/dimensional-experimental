-- | This module provides automatic differentiation for Quantities.

{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module Numeric.Units.Dimensional.AD (FAD, diff, Lift (lift), undim, todim) where

import Numeric.Units.Dimensional (Dimensional (Dimensional), Quantity, type (/))
import Numeric.Units.Dimensional.Coercion
import Numeric.AD (AD, Mode, auto, Scalar)
import qualified Numeric.AD (diff)
import Numeric.AD.Mode.Forward (Forward)

-- | Unwrap a Dimensional's numeric representation.
undim :: Quantity d a -> a
undim = coerce
todim :: a -> Quantity d a
todim = coerce

type FAD tag a = AD tag (Forward a)

-- | @diff f x@ computes the derivative of the function @f(x)@ for the
-- given value of @x@.
diff :: Num a
     => (forall tag. Quantity d1 (FAD tag a) -> Quantity d2 (FAD tag a))
     -> Quantity d1 a -> Quantity ((/) d2 d1) a
diff f = todim . Numeric.AD.diff (undim . f . todim) . undim

-- | Class to provide 'Numeric.AD.lift'ing of constant data structures
-- (data structures with numeric constants used in a differentiated
-- function).
class Lift w where
  -- | Embed a constant data structure.
  lift :: (Mode t) => w (Scalar t) -> w t

instance Lift (Quantity d)
  where lift = todim . auto . undim
