-- | This module provides automatic differentiation for Quantities.

{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE Rank2Types #-}

module Numeric.Units.Dimensional.AD (diff, Lift (lift)) where

import Numeric.Units.Dimensional (Dimensional (Dimensional), Quantity, Div)
import Numeric.AD (AD, Mode, auto, Scalar)
import Numeric.AD.Mode.Forward (Forward)
import qualified Numeric.AD (diff)

-- | Unwrap a Dimensional's numeric representation.
undim :: Dimensional v d a -> a
undim (Dimensional a) = a

-- | @diff f x@ computes the derivative of the function @f(x)@ for the
-- given value of @x@.
diff :: (Num a, Div d2 d1 d3)
     => (forall tag. Quantity d1 (AD tag (Forward a)) -> Quantity d2 (AD tag (Forward a)))
     -> Quantity d1 a -> Quantity d3 a
diff f = Dimensional . Numeric.AD.diff (undim . f . Dimensional) . undim

-- | Class to provide 'Numeric.AD.lift'ing of constant data structures
-- (data structures with numeric constants used in a differentiated
-- function).
class Lift w where
  -- | Embed a constant data structure.
  lift :: (Mode t) => w (Scalar t) -> w t

instance Lift (Dimensional v d)
  where lift = Dimensional . auto . undim
