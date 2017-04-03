This repository contains the code for my Bachelor thesis. There is some (incomplete) documentation of C++ code in the `docs` directory and `README.pdf` contains some explanation of the underlying physical theory. 

The following list contains some questions I could find no answer to yet. I'd really appreciate it if you could have a look at them and maybe help me out.

Questions
---------

1) In RPA, we use perturbation potential of the form `V*exp(-i * omega * t + eta * t)`. Why is `V` taken to be diagonal in position representation? This assumtion seems logical for some cases. Consider, for example, constant external electric field E. Then `V(r) = -q * E . r`, which is, by definition, diagonal. The same holds for applied magnetic field `B(r)`. But does this hold in general?

__Answer__: I spoke to Mikhail Katsnelson today, he confirmed that `V` is usually diagonal.

2) Have a look at the derivation in footnote 3. What is the Fourier transform of it? Factor exp(eta * t) results in divergencies in Fourier integral... Or is there a mistake in my calculations? Maple, for example is able to calculate the transform of `exp(-i * omega * t) * (1 + eta * t)`. Here, I used the Taylor expansion of `exp`. The result is `2 * Pi * (I * eta * Dirac(1, omega) + Dirac(omega))`. This is a really nice, because then we can neglect the first delta function under the assumption that `eta` is small. Unfortunately, I failed to derive this expression myself and I don't really trust my knowledge of Maple in this case.

__Answer__: There is a rather nice description of this concept in "Causality and Dispertion Relations" book by H.M. Nussenzveig. The method's description should now be correct -- see updated PDF.

3) ...

