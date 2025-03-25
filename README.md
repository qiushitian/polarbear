# Polarbear, the simplistic reduction scripts for ARCTIC on the APO 3.5m telescope

*"How boring would the arctic be if there were no polarbears!"* —Chris

Author: Qiushi (Chris) Tian

Last updated: 2023-03-25

## Dependencies

- NumPy
- Astropy
- ccdproc

## What does it do?

Currently, there are only three functions:

- `trim_overscan`. It only *trims* the overscan areas but does not do overscan subtraction.

- `make_master_bias`. It creates a sigma-clipped mean master bias.
There are parameters you can tweak in the `ccdproc.combine` call.

- `flatfield`. It performs flat-correction to science frames. Requires a master bias frame.
*Note: `flatfield` will call `trim_overscan` to trim overscan after flatfielding.
The idea is to make make the master bias* without *trimming overscan,
so that you don't need to trim the science frames' overscan one by one.*

There is also a `if __name__ == '__main__'` section.
So, if you put lines of code there and run the file, it can run like a script.
There are some examples there.
