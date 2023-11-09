## handyFunctionObjects

### Installation
Copy the directory somewhere, ensure the correct version of OpenFOAM is sourced, and run wmake

``` bash
git clone https://github.com/Ben-Malin/handyFunctionObjects.git
cd handyFunctionObjects && wmake
```

This will compile a new library, libhandyFuncs.so, in your $FOAM_USER_LIBBIN directory.

To make use of the function objects, this library will need to be included in your controlDict, or passed to the postProcess application using the -lib flag

In controlDict:
``` c++
libs
(
    "libhandyFuncs.so"
);
```

For postProcessing (note that objects which make use of a turbulence model require postProcess to be called through a solver):
```bash
simpleFoam -postProcess -func SSTBlending -latestTime -lib libhandyFuncs.so
```

### Function Objects
#### SSTBlending

The kOmegaSST turbulence model uses a blending function, to apply the $k-\omega$ model in the near wall region, and the $k-\epsilon$ model in the freestream.  
Sometimes, it might be interesting or useful to actually know where this blending takes place.
This function object simply writes out the F1 blending function that is used to select the model coefficients.

Once installed, the functionObject can be run after the simulation using the above postProcess command, or by including the following in the functions section of controlDict:

```c++
functions
{
    sstBlending
    {
        type SSTBlending;
        libs ("libhandyFuncs.so")
        executeControl writeTime;
        writeControl   writeTime;
    }
}
```
