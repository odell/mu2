# mu2

`mu2` is a set of two-body scattering and bound-state codes for arbitrary
interactions (except Coulomb).

`System` is the primary class that the user interacts with. Most importantly, it
takes an `Interaction` instance as input --- the interaction that models the
physics between the two particles.

`Interaction` is a combination of at least 2 pieces: the "long-range potential",
defined in coordinate space, and the `Counterterm`. The long-range potential is
regulated at the short distance, `R`.

`Counterterm` defines the leading-order and next-to-leading-order "contact"
interactions that allow the user to tune the short-distance physics.

The supported regulation schemes are nonlocal, semi-local, and local.
