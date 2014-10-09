{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}

-- | Provides instances of @Data.AEq.AEq@ from the ieee754 library.
module Numeric.Units.Dimensional.AEq where

import Numeric.Units.Dimensional (Dimensional (Dimensional), Quantity)
import Numeric.Units.Dimensional.LinearAlgebra.Vector (Vec (ListVec))
import Data.AEq


deriving instance AEq a => AEq (Quantity d a)

instance (Floating a, AEq a) => AEq (Vec ds a)  -- CPos et al
  where
    -- ListVec xs === ListVec ys = and $ zipWith (===) xs ys
    ListVec xs ~== ListVec ys = and $ zipWith (~==) xs ys
