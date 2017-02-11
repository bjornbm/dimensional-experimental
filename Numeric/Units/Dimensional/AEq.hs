{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

-- | Provides instances of @Data.AEq.AEq@ from the ieee754 library.
module Numeric.Units.Dimensional.AEq where

import Numeric.Units.Dimensional (Dimensional (Dimensional), Quantity)
import Numeric.Units.Dimensional.Coercion
import Numeric.Units.Dimensional.LinearAlgebra.Vector (Vec (ListVec))
import Numeric.Units.Dimensional.LinearAlgebra.Matrix (Mat (ListMat))
import Numeric.Units.Dimensional.DK.AEq
import Data.AEq


instance AEq a => AEq (Vec d n a)
  where
    ListVec xs ~== ListVec ys = xs ~== ys

instance AEq a => AEq (Mat d r c a)
  where
    ListMat xs ~== ListMat ys = xs ~== ys
