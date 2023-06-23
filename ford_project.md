---
project: ALPS
project_github: https://github.com/danielver02/ALPS
author: Kristopher Klein, Daniel Verscharen
email: kgklein@arizona.edu, d.verscharen@ucl.ac.uk
srd_dir: src
page_dir: docs
output_dir: docs
fpp_extensions: f90
base_url: http://www.alps.space
ordered_subpage: INSTALL.md
                 tutorial.md
                 input.md
                 output.md
display: public
         protected
         private
source: false
graph: true
incl_src: false
src_dir:  ./src
          ./interpolation
          ./distribution

---

<img src="./sub_doc/ALPS_logo.png" alt="drawing" width="200"/>
<img src="./sub_doc/qrcode_alps_github.png" alt="drawing" width="200"/>

ALPS is a parallelised numerical code that solves the Vlasov-Maxwell dispersion
relation in hot (even relativistic) magnetised plasma. ALPS allows for any
number of particle species with arbitrary gyrotropic equilibrium distribution
functions supporting waves with any direction of propagation with respect to
the background magnetic field.

If you use the code for a science publication, please provide the code website
[github.com/danielver02/ALPS](https://github.com/danielver02/ALPS) in the acknowledgements of your publication and cite the code paper:

[Verscharen, D., Klein, K. G., Chandran, B. D. G., Stevens, M. L., Salem, C. S.,
and Bale, S. D.: ALPS: the Arbitrary Linear Plasma Solver, J. Plasma Phys. 84,
905840403, 2018, doi: 10.1017/S0022377818000739](http://doi.org/10.1017/S0022377818000739)
