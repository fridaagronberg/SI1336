# Simulation polymers in 3D using random walks and freely joint chains

## Problem description

From Canvas: <i> Polymers: compare a 3D random walk on a grid with a freely jointed
chain (chain with fixed-length links with random orientation), both
non self-avoiding as well as self-avoiding (note that you need to use
spheres with a chosen radius for the self-avoiding freely jointed
chain). Take care to generate a uniform distribution of angles for the
freely-jointed chain.</i>

## Parameters

**nsteps** - number of steps: 10-20000
**self_avoiding** - on/off
**can walk backwards**
**radius** of the sphere that checks if the freely joint chain is avoiding: 0.001-1
**parameters in random_int_generator**: primarily how does it affect random walk,
is it possible to use the same function to decide angle of the freely joint chain?

## "Measurable" quantities
**distance** - from origo until endpoint, rms over multiple simulations with the same parameters
- discard or use walks that cross themselves, when calculating rms? both?
- for spec nsteps plot number of walks that are a certain distance from origo when stopping (either crashing in to itself or finished)

**variance and fluctuation estimate**

**How uniform the walks and chains are**
- depending on random number generator
- nsteps

## Discretization

From canvas: *comparison of relevant methods/parameters, motivated choices*

Fulfilled by:
- discussion of the parameters effect and the choice of random_int_generator
- maybe, statistical test of random_int_generator randomness using \chi^2
- no time discretization in this case but should be fine

### Small comments
**self-avoiding** = cannot cross itself
Odds that it returns to origo?

## Models

### RandomWalk

### FreelyJointedChain

Även om längden är "1" så måste avståndet mellan punkterna vara mindre (finer grid) annars kan man bara ta 14steg (kub)
Borde vara såpass fint att om den tar ett rakt eller ett diagonalt steg inte påverkar totala längden så mycket.

Känns fel att modellera som diskreta steg i ett grid... jag trollar... är bara att räkna vektoriellt...

Modellera med sfäriska koordinater och generera random phi och theta, undersöka random_int_generator ; ej på varandra beroende intergers.
