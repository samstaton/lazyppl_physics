{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TemplateHaskell       #-}


module Physics where
import Apecs.Physics
import Apecs.Physics.Gloss 

import Graphics.Gloss.Export (exportPicturesToGif,GifLooping(LoopingForever))

import Control.Monad (replicateM, when) 

import LazyPPL
import LazyPPL.Distributions (uniformbounded, normalPdf)

import Numeric.Log(Log(Exp),ln)
import System.Random (setStdGen, mkStdGen)
import Graphics.Matplotlib hiding (density)

import System.IO.Unsafe (unsafePerformIO)
import Data.List (findIndex)
import Data.Monoid (Product(Product))

{-- A simple 2d physics problem based on the Chipmunk / Apecs library
    Inspired by a tutorial example in Anglican. See e.g.:
    https://github.com/hongseok-yang/probprog19/blob/master/Lectures/Lecture4/PhysicsLec4.clj
    Thanks to Alexander Bai for adapting this from the Haskell example in Sam's OPLSS course to LazyPPL.
 --}

{-- Set up Apecs Physics --} 
makeWorld "World" [''Physics, ''Camera]

{-- Initialize the world a ball, a cup, and bumpers of given position and angle.
    Return the ball. --} 
initialize bumpers = do
  set global ( Camera (V2 0 1) 60
             , earthGravity )

  -- The bumpers 
  mapM (\(x,y,theta) -> do
          lineBody <- newEntity (StaticBody, Angle (theta), Position (V2 x y))
          newEntity (Shape lineBody (hLine 2), Elasticity 0.3)
       ) bumpers

  -- The cup
  lineBody <- newEntity (StaticBody, Position (V2 4.5 (-5)))
  newEntity (Shape lineBody (hLine 1), Elasticity 0.1)  
  lineBody <- newEntity (StaticBody, Position (V2 4 (-4.5)))
  newEntity (Shape lineBody (vLine 1), Elasticity 0.1)  
  lineBody <- newEntity (StaticBody, Position (V2 5 (-4.5)))
  newEntity (Shape lineBody (vLine 1), Elasticity 0.1)  

  -- The ball
  ball <- newEntity (DynamicBody, Position (V2 (-4.5) 5))
  newEntity (Shape ball (cCircle 0.2), Density 1, Elasticity 0.9)
  return ball

-- Criterion for ending the simulation:
-- the ball coordinates falls below the bottom of the window
-- or into the cup
endCriterion :: (Double,Double) -> Bool
endCriterion (x,y) = y< -5 || (y < -4.5 && x > 4 && x < 5)


-- Make an animated gif from a sequence of pictures using Gloss
exportGif filename pics = do
   putStrLn $ "Plotting GIF to " ++ (show filename) ++ "..."
   exportPicturesToGif 2 LoopingForever (400,400) (makeColor 1 1 1 0) filename (\t -> pics !! (floor t)) [0..(fromIntegral $ length pics - 1)]
   putStrLn $ "Done."

-- Plot the given trajectories to an SVG file
plotCoords :: String -> [[(Double,Double)]] -> Double -> Double -> Double -> IO ()
plotCoords filename  xyss ymin ymax alpha = 
    do  putStrLn $ "Plotting " ++ filename ++ "..."
        file filename $ mp # "plot.axis('off')" % mp # "colors = plot.cm.gnuplot2(np.linspace(0,1," # (floor ((fromIntegral $ length xyss) * 1.5) :: Integer) # ")) " % foldl (\a i -> a % plot (map fst (xyss !! i)) (map snd (xyss !! i)) @@ [o1 "o-", o2 "color" (lit ("colors[" ++ (show (length xyss - i)) ++ "]")), o2 "linewidth" (1.0 :: Double), o2 "alpha" alpha, o2 "ms" (0 :: Int)]) mp [0..(length xyss -1)] % (plot ([4,4,5,5] :: [Double]) ([-4,-5,-5,-4] :: [Double]) @@ [o1 "o-", o2 "c" "black", o2 "ms" (0 :: Int)] % xlim (-5 :: Int) (5.1 :: Double) % ylim ymin ymax) 
        putStrLn "Done."
        return ()


-- Given the bumper positions,
-- return the trajectory of the ball
-- (a list of (x,y) coordinates)
getTrajectory :: [(Double,Double,Double)] -> IO [(Double,Double)]
getTrajectory bumpers =
  do w <- initWorld
     runWith w $ do
       b <- initialize bumpers
       xys <- run b 500 
       return xys
  where
    -- run produces the trajectory until the endCriterion
    -- or until the fuel runs out
    -- (in case the ball somehow gets stuck)
    run b fuel = do
      stepPhysics (1/60)
      (Position (V2 x y)) <- get b
      if fuel < 0 then return [(x,y)] else
        if endCriterion (x,y) then return [(x,y)] else
        do
           rest <- run b (fuel - 1)
           return $ (x,y) : rest
   
-- Given the bumper positions,
-- return pictures showing the scenario over time.
getPics :: [(Double,Double,Double)] -> IO [Picture]
getPics bumpers =
  do w <- initWorld
     runWith w $ do
       b <- initialize bumpers
       xypics <- mapM (\_ -> do
                       stepPhysics (1/60)
                       (Position (V2 x y)) <- get b
                       pic <- draw 
                       return (x,y,pic)
                   ) [1..1000]
       let (Just n) = findIndex (\(x,y,_) -> endCriterion (x,y)) xypics
       return $ map (\(_,_,pic) -> pic) $ take (n + 50) xypics
  where
    draw = do
      -- draw the current scene
      pic <- foldDrawM drawBody
      let cam = Camera (V2 (5) (-5)) 20
      return $ cameraTransform cam pic 

-- The model: Pick some bumper positions and angles uniformly
-- Run the 2d physics
-- Then observe that the ball lands more-or-less in the cup.
-- Return a posterior distribution over bumper positions and angles.
model :: Meas [(Double,Double,Double)]
model = do
  -- Pick the positions and angles of two bumpers
  -- Prior is uninformative, uniformly distributed
  bumpers <- sample $ replicateM 2 $ do
                      x <- uniformbounded (-5) 5
                      y <- uniformbounded (-5) 5
                      theta <- uniformbounded 0 pi
                      return (x,y,theta)
  -- Run the physics simulation.
  -- The last point is where the ball either leaves the scene
  -- or lands in the cup. 
  let (x,y) = unsafePerformIO $
              do { pos <- getTrajectory bumpers ; return $ last pos }
  -- Observe the ball is roughly in the cup at the end of the trajectory
  score $ normalPdf 4.5 0.2 x 
  -- Return the bumper positions
  return bumpers

-- Run the model with MH
runModelMH n trajfile anifile = do
  -- Use MH with p=0.15 (0.15 ~= 1/6, and approx six samples needed)
  samples <- mh 0.15 model
  -- The ball enters the cup when the likelihood is > 0.0876
  let p = Product $ Exp $ ln $ 0.0876415
  let (Just i) = findIndex (\(_,w) -> w > p) samples
  putStrLn $ "[MH] Ball entered the cup after " ++ show i ++ " samples."
  -- Output illustrations to files:
  when (trajfile /= "") $ do
     -- Extract the trajectories of the first n samples
     xyss <- mapM getTrajectory (map fst $ take n samples)
     -- Plot the trajectories together
     plotCoords trajfile xyss (-5) 5 0.01
  when (anifile /= "") $ do
    -- Plot the last trajectory
    pics <- getPics $ fst $ samples !! n
    exportGif anifile pics

-- Run the model with brute force search
runModelBruteForce n trajfile = do
  samples <- weightedsamples model
  -- The ball enters the cup when the likelihood is > 0.0876
  let p = Exp $ ln $ 0.0876415
  let (Just i) = findIndex (\(_,w) -> w > p) samples
  putStrLn $ "[BF] Ball entered the cup after " ++ show i ++ " samples."
  -- Extract the trajectories of the first n samples
  xyss <- mapM getTrajectory (map fst $ take n samples)
  -- Plot the trajectories together
  plotCoords trajfile xyss (-5) 5 0.1


main :: IO ()
main = do
     -- The distribution is multimodal,
     -- so let's restart four times for four examples.
     -- (Statistically better to use mhirreducible.)
     -- We plot the trajectories of the last run.
     setStdGen (mkStdGen 1)
     runModelMH 1000 "" "physics1.gif"
     setStdGen (mkStdGen 42)
     runModelMH 1000 "" "physics2.gif"
     setStdGen (mkStdGen 47)
     runModelMH 1000 "" "physics3.gif"
     setStdGen (mkStdGen 57)
     runModelMH 1000 "physics4.svg" "physics4.gif"
     -- We also plot what brute force looks like with 4000 trials
     runModelBruteForce 4000 "physics5.svg"