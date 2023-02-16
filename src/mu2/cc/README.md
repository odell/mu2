# mu2 C Implementations

`mu2` offers the capability to use C code for scattering calculations. The
relevant files are contained here in the `cc` directory.

In order to call the appropriate methods of the System instance, denoted with a
trailing `_fast`, the user needs to compile the `libkcd.so` shared library. The
`makefile` exists to make this process easier, via

```
make libs
```

The shared library is placed in the `$HOME` directory for convenience.
