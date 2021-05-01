{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NoStarIsType #-}

module Numeric.Units.Dimensional.LinearAlgebra.VectorAD where

import qualified Prelude
import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional (Dimensional (Dimensional))
import Numeric.Units.Dimensional.LinearAlgebra.Vector (Vec (ListVec), elemAdd, scaleVec)
import Numeric.AD (auto, diffF')
import Numeric.Units.Dimensional.AD


-- | If @f@ is a function of a quantity that returns a 'Vector', then
-- @diff f@ is a function of the same type of quantity that returns
-- the first derivative of the result.
diffV :: (Num a, d' ~ (/) d t)
      => (forall tag b . b ~ FAD tag a => Quantity t b -> Vec d n b)
      -> Quantity t a -> Vec d' n a
diffV f = snd . diffV' f


unfzip as = (fmap fst as, fmap snd as)

-- | Like 'diffV' but returns a pair of the result and its first derivative.
diffV' :: (Num a, d' ~ (/) d t)
       => (forall tag b . b ~ FAD tag a => Quantity t b -> Vec d n b)
       -> Quantity t a -> (Vec d n a, Vec d' n a)
diffV' f x = (ListVec ys, ListVec ys')
  where
    (ys,ys') = unfzip $ diffF' (unvec . f . todim) (undim x)
    unvec (ListVec xs) = xs


--primalV :: Num a => Vec ds (AD tag a) -> Vec ds a
--primalV (ListVec xs) = ListVec (fprimal xs)


-- Linearizing
-- -----------

-- TODO
-- -- | @applyLinear@ converts a pair of a vector and its derivative w r t a
-- --   variable (e g time) into a function  linearized about the original vector
-- --   at @t=0@. Then the function (which should be independent of the variable,
-- --   but see 'applyLinearAt') is evaluated and the "new" vector/derivative
-- --   pair is reconstructed from the result.
applyLinear :: forall d1 d1' d2 d2' n m a t
            . (d1 ~ (*) t d1', d2' ~ (/) d2 t, (/) d1 d1' ~ t, Num a)
              --  HZipWith DivD ds ds' ts, Homo ts t  -- Necessary to infer t (the dimension w r t which we are differentiating).
            => (forall tag b . b ~ FAD tag a => Vec d1 n b -> Vec d2 m b)
            -> (Vec d1 n a, Vec d1' n a) -> (Vec d2 m a, Vec d2' m a)
applyLinear f (v,v') = applyLinearAt (const f) (_0 :: Quantity t a) (v,v')

-- | 'applyLinearAt is analogous to 'applyLinear' but should be used when
--   the function is also dependent on the variable w r t which the vector
--   is linearized.
applyLinearAt :: (d1 ~ (*) t d1', d2' ~ (/) d2 t, Num a)
              => (forall tag b . b ~ FAD tag a => Quantity t b -> Vec d1 n b -> Vec d2 m b)
              -> Quantity t a -> (Vec d1 n a, Vec d1' n a) -> (Vec d2 m a, Vec d2' m a)
applyLinearAt f t (v,v') = diffV' (\t' -> f t' (lift v `elemAdd` scaleVec (t' - lift t) (lift v'))) t


-- Lifting
-- -------

-- | Lift the elements of a vector to 'AD.AD's.
instance Lift (Vec d n) where lift (ListVec xs) = ListVec (map auto xs)
