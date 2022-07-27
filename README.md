## 2d physics simulation in LazyPPL

Attempt to add physics simulation as an example of lazyppl

See the original file at Sam Staton's OPLSS course [[link]](https://www.cs.uoregon.edu/research/summerschool/summer19/topics.php#Staton)

See LazyPPL [here](https://lazyppl.bitbucket.io/)

## Installation

The system uses Haskell stack.
You need to [install stack](https://docs.haskellstack.org/en/v1.1.2/install_and_upgrade/) first if you want to use the system in the
standard way.

To build, type
``stack build``.
This may take some time (>1 hour) if it is your first time ever using stack.

To run, type
``stack run`` or ``stack run physics-exe``

## Note

The original OPLSS sample code fails to compile in macOS (at least my labtop).  A similar problem is seen [here](https://github.com/bravit/hid-examples/issues/7). LazyPPL switched to use [matplotlib](https://hackage.haskell.org/package/matplotlib) in their development. A similar problem is also seen in other haskel probablistic programming tutorials ([like this one](https://github.com/ccshan/prob-school)).
