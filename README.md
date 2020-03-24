# Smallptsampling

We use the [smallpt](https://www.kevinbeason.com/smallpt/) source code
to show off different sampling techniques. We use the original smallpt
as a baseline, and then modify the sampler to be more and more complex.


- `baseline.cpp` is the code from smallpt, `ClangFormat`ted.
- `baseline-serial.cpp` is the code from smallpt with parallelism removed.
  This will help to try and use other samplers.
- `baseline-traced.cpp` is the code from smallpt, serialized, then built
  using the `Trace` concept. This is to explore whether using `Trace` allows
  us to get better samples.
