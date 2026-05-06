# QuiverCombinatoricsTools
 
AGQ computing project by Mia, Tudor, and Emanuel. See QuiverCombinatorics.pdf for our aims and goals.

We wish to extend the SAGEmath package [QuiverTools](https://github.com/QuiverTools/QuiverTools) written by Pieter Bielmans, Hans Franzen, and Gianni Petrella. This is the extension of the package. You can read the documentation of QuiverTools as

* [a webpage](https://sage.quiver.tools)
* [a pdf](https://sage.quiver.tools/documentation.pdf)

A more detailed user guide is in the works.

# Instructions

You can install it by going into your sage environment, and ensuring you can get to this repository. Because it is private, you will have to install an ssh key! Make sure an .ssh/authorized_keys directory exists, then run

``nano ~/.ssh/authorized_keys``

and add the following line

``ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIKHt2WDoUV+R+t2sjWwsRrhYqbJiZJJsoagfGoNnlDXb QuiverCombinatoricsTools``

then exit nano. Now, enter the environment with your sage install, and run pip

``pip install git+https://github.com/emanuel-roth/QuiverCombinatoricsTools``

in any sage code, use `from quivercombinatorics import *` to get started.