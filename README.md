# QuiverCombinatoricsTools
 
AGQ computing project by Mia, Tudor, and Emanuel.

We hope to extend the sage package [QuiverTools](https://github.com/QuiverTools/QuiverTools) written by Pieter Bielmans, Hans Franzen, and Gianni Petrella. We do this by copying the package here and extending this repository.

You can read the documentation of QuiverTools as

* [a webpage](https://sage.quiver.tools)
* [a pdf](https://sage.quiver.tools/documentation.pdf)

# Instructions

Since this is a private repository for now, you can access it as follows. In your unix bash terminal, add an .ssh folder if you don't have it already

``mkdir -p ~/.ssh``

Then open nano to make authorized_keys

``nano ~/.ssh/authorized_keys``

and add this key to the file

``ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIKHt2WDoUV+R+t2sjWwsRrhYqbJiZJJsoagfGoNnlDXb QuiverCombinatoricsTools``

Then open nano in the ssh config

``nano ~/.ssh/config``

and add this

```
Host github      
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_ed25519
  IdentitiesOnly yes
  ```

Where you have sage installed, run 

``pip install git+ssh://git@github.com/emanuel-roth/QuiverCombinatoricsTools.git``

Then, in any sage code, use `from quivercombinatorics import *` to get started.