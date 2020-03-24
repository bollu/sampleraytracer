# Smallptsampling

We use the [smallpt](https://www.kevinbeason.com/smallpt/) source code
to show off different sampling techniques. We use the original smallpt
as a baseline, and then modify the sampler to be more and more complex.


- `baseline.cpp` is the code from smallpt, `ClangFormat`ted. This code
  uses the same randomness across each `y` axis.
- `baseline-serial.cpp` is the code from smallpt with parallelism removed.
  This will help to try and use other samplers.
- `baseline-traced.cpp` is the code from smallpt, serialized, then built
  using the `Trace` concept. This is to explore whether using `Trace` allows
  us to get better samples.
- `baseline-same-randomness.cpp` is the code that uses the _same_ randomness
  for each pixel. This causes us to see "streaks".

- `baseline-xy-randomness.cpp` is the code that uses the different randomness
  for each pixel (in x and y axis).


# Story

first made `baseline-serial` to check if parallelism has nondeterminism
or if it was equivalent to serial semantics.

Then checked effect of type of randomness by creating baseline-same-randomness
and baseline-xy-randomness.
