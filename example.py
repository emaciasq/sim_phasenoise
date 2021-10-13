import os
try:
    # Needed if you want to produce plots in simobserve_phase
    import vis_tools
    execfile('{}/CASA_vis_tools.py'.format(os.path.dirname(vis_tools.__file__)))
except:
    print('Cannot find vis_tools in your path. The module vis_tools and ' +
    'accompanying script CASA_vis_tools.py is needed to run some of the ' +
    'plots. You can find them at https://github.com/emaciasq/vis_tools')

execfile('/path/to/sim_phasenoise/sim_phasenoise.py')

fconf = '/path/to/array/conf/X.cfg'     # Configuration file
pwv = 0.5                               # pwv (mm)
obstime = str(3.0 * 3600.0)+'s'         # Total observing time
integration = '8s'                      # Integration time of each data point
prms_0 = ['3deg', '8deg', '12deg']      # Phase rms
prms_timescale = ['24s', '72s', '120s'] # Timescales to introduce phase noise
project = 'your_simulation'             # Name for project
model_image = 'model.fits'              # Path to model fits image
bandwidth = 8.0                         # Bandwidth (GHz)

simobserve_phase(fconf, pwv, obstime, model_image = model_image,
                 bandwidth = bandwidth, prms_0 = prms_0,
                 integration = integration, prms_timescale = prms_timescale,
                 plots = True, project = project)

# The simulated measurement set with only thermal noise will be at:
ms_thermalnoise = project+'/'+project+'.'+os.path.basename(fconf)[:-3]+'noisy.ms'

# The simulated measurement set with phase noise will be at:
ms_phasenoise = project+'/'+project+'.'+os.path.basename(fconf)[:-3]+'noisy_phasenoisy.ms'


