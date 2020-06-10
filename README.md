Whack-A-Mole

## Background Story
The character is a farmer who is devastated by mole pestilence in their farm. Now the farmer wants to start a new seed to plant a new plant in the garden, but first the farmer should remove all the moles to guarantee success of the crops this season !!!


## Gameplay
* There will be several holes (depending on how many marker)
* The mole will appear randomly out of the hole
  * Many mole at many hole
* How long is the game?
  * to be determined (2~3 min?)
* How to calculate the score
  * predetermined number: to be determined (-2?)
  * penalty, cannot hit the appeared mole to be determined (-1?), 
  * if we have time → (-1 (easy), -2(med), -3(hard)) 



## Mechanics
(How to play the game)
Hitting mole mechanics, Is it 
* on-screen mouse click mechanics (on pre-recorded image/ video mode)
* ~~marker collision (the hand also has marker)~~
* marker loss (on camera mode)

1. Firstly, the game has to detect how many (e.g six) markers were set up before starting the game? 
2. Push starting button to start (exercise image)
3. Do the mole appear at the same time? Or one in a time, how long 
   1. random going up time
   2. Interval of mole appear up-top 
   3. (to be decided) seconds
4. Start again button

## Assets Needed
* Start/Restart button (essential feature) simple 2D
* Hole (essential feature) (3D)
* Mole (essential feature) (3D)
* Land/Grass (background) (additional feature)


## Team Member
* Fariz Ikhwantri
* Chen Pinjung
* Akira Watanabe
* Chaijirawiwat Chawit
* Kaneko Mayu

# Schedule
6/11 progress update
6/18 presentation

## Work Assignment
1. asset rendering [1.1] (Akira)
   1. mole
   2. hole (imagined base object: half Torus)
2. marker detection and table (image and video) [1.3] (Chawit)
   1. hole_list = [hole1 hole2 .. hole6] ← (ID, perspective stored here) 
3. game mechanics
   1. state saving template [3] (Fariz): state saving (how many moles at frame i-th currently appear)?  (knowing what previous code has done)
      1. HoleState (repurpose/use/include MarkerTracker)
      2. MoleState
      3. instead of integrating all of game mechanics implementation later, each assigned member use provided hole and mole state objects to modified the state
   2. transition animation [1.2] (Kaneko)
      1. pop-up (hole and mole state)
      2. pop-down (hard because there will be no marker) (hole and mole state)
      3. hit???? (hard because there will be no marker) (mole state and then hole state)
   3. marker redetected again
   4. marker loss detection? [1.2 → 2.1]  bool ifMarkerExist[6]
   5. random mole appearance? [2.2] (Chen) holesState[i] = [0 0 1 0 1 0] 
   6. scoring mechanics [4]
4. integration (should we do this together ?? because it will need at ) if we are doing call from provided 


(Chawit) Wed 10/6 → (Fariz) Sat 13/6→ (Chen) 15/6→  integration

(Akira) Wed 10/6 → (Kaneko) → integration
(Fariz) Sat 13/6    ↗ 

Presentation 






