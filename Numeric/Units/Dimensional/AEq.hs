{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}

-- | Provides instances of @Data.AEq.AEq@ from the ieee754 library.
module Numeric.Units.Dimensional.AEq where

import Numeric.Units.Dimensional (Dimensional (Dimensional), Quantity)
import Numeric.Units.Dimensional.LinearAlgebra.Vector (Vec (ListVec))
import Numeric.Units.Dimensional.LinearAlgebra.Matrix (Mat (ListMat))
import Data.AEq


deriving instance AEq a => AEq (Quantity d a)

instance (Floating a, AEq a) => AEq (Vec ds a)
  where
    ListVec xs ~== ListVec ys = xs ~== ys

instance (Floating a, AEq a) => AEq (Mat ds a)
  where
    ListMat xs ~== ListMat ys = xs ~== ys
