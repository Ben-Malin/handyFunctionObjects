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
simpleFoam -postProcess -latestTime -lib libhandyFuncs.so
```

### Function Objects
#### SSTBlending

The $k\omega SST$ turbulence model uses a blending function, to apply the $k-\omega$ model in the near wall region, and the $k-\epsilon$ model in the freestream.  
Sometimes, it might be interesting or useful to actually know where this blending takes place.
This function object simply writes out the F1 blending function that is used to select the model coefficients.

Once installed, add the following to the functions section of your controlDict:

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

And the blending field will be saved during your simulation.
If you want to output the field for a simulation you've already run, you can use the simpleFoam -postProcess command as above
