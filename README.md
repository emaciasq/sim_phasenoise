# sim_phasenoise
Tools to include atmospheric noise in simulations of radio-interferometric
observations using CASA.

`sim_phasenoise` extends the capabilities of CASA's
`simobserve` task, which currently can only include thermal noise. The phase
noise is introduced using a method developed by Ed Fomalont: a "de-calibration"
table is populated with phase noise, which is then  applied to the data using
`applycal`.

The script can currently include phase noise at multiple timescales, although in
a *step-wise* way that is not fully realistic. The phase noise is also
introduced independently for each antenna (so the noise for nearby antennas will
not be correlated). Despite these current limitations, `sim_phasenoise` can
provide more realistic simulations than the ones produced simply by
`simobserve`.

## Basic usage

This script should be run inside CASA. The firt step is to load the scripts:

```python
CASA <X>: execfile('/your/path/to/sim_phasenoise.py')
```

The most direct way of applying phase noise to a simulation would be using the
`apply_phasenoise` function. This function takes as inputs a CASA measurement
set (not necessarily from a simulation, in fact) and an array configuration
file, and creates a new measurement set with phase noise included:

```python
CASA <X>: apply_phasenoise('/path/to/measurement/set', '/path/to/array/configuration/file')
```

The function will create a mock gain calibration table called atmos.cal, which
will be populated with phase noise and then applied to the data.

## Atmospheric phase rms

The phase noise is computed based on the phase rms produced by the spatial
structure function measured during the ALMA Long Baseline campaign between 2012
and 2014 by Matsushita et al. (2017, PASP, 129, 035004). This spatial structure
function uses two power-laws, with a slope of 0.8 up to baselines of 10 km, and
0.25 beyond 10 km. By default, it is scaled to give a phase rms of 30 degrees at
10 km, which is a bit better than the nominal go-nogo conditions at ALMA (1
radian). However, note that phase referencing is always used to calibrate the
phases, thus reducing their rms further. The parameter `prms_0` can be used to
set a lower (or higher) phase rms at 10 km. `apply_phasenoise` has other
optional input parameters that allow, among other things, to apply different
phase rms at different timescales. See the docstring for more information.

## simobserve_phase

`sim_phasenoise` also has a wrapper function called `simobserve_phase`, which
runs the CASA task `simobserve` and then the function `apply_phasenoise`.
`simobserve_phase` takes as input parameters the array configuration file, the
precipitable water vapor (`pwv`) to be used in `simobserve`, the total observing
time (`obstime`), and the input model. The model can be provided using a fits
image:

```python
CASA <X>: simobserve_phase('/path/to/array/configuration/file', pwv, obstime,
                           model_image = '/path/to/model/image')
```

The model can also be set by providing a list of Gaussian components. See the
docstring for more info.

`simobserve_phase` can produce a few plots if `plots` is set to True. In order
to create these plots you will need to have in your path the vis_tools module
(https://github.com/emaciasq/vis_tools), and will need to run:

```python
CASA <X>: import vis_tools
CASA <X>: execfile('{}/CASA_vis_tools.py'.format(
              os.path.dirname(vis_tools.__file__)))
```

## Example

The script `example.py` shows an example of how to run `simobserve_phase`.
