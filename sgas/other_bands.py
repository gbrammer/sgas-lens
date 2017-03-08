def prep():
    """
    Rerun the WFC3 calibration pipeline to flatten the (potentially)
    variable ramps
    """
    import matplotlib as mpl
    mpl.rcParams['backend'] = 'agg'

    import glob
    import os
    
    import stsci.tools
    from sgas import reprocess_wfc3
    
    files=glob.glob('*raw.fits')
    reprocess_wfc3.show_ramps_parallel(files, cpu_count=4)

    files=glob.glob('*raw.fits')
    reprocess_wfc3.reprocess_parallel(files)
    
    