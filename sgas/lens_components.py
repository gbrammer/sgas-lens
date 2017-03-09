#!/usr/bin/env python
# encoding: utf-8
"""
Code for fitting individual components of the lens+source system

"""

# The ID of the *lens* in my segmentation file, could be different in
# yours
LENS_SEG_ID = 253

def extract_beams():
    """
    Get beam cutouts from the grism exposures
    """
    import glob
    import drizzlepac # to avoid memory issues with GroupFLT...
    import numpy as np
    
    import grizli
    import grizli.prep
    
    # Find the grism exposures
    files=glob.glob('*_flt.fits')
    info = grizli.utils.get_flt_info(files)
    grism_files = list(info['FILE'][info['FILTER'] == 'G141'])
    np.save('all_grism_files.npy', [grism_files])
    
    # Initialize GroupFLT object
    # This should read GrismFLT files and static contam model
    # you generated in grism.py
    grp = grizli.multifit.GroupFLT(grism_files=grism_files,
                                   direct_files=[], 
                                 ref_file='sdssj0851+3331-f160w_drz_sci.fits',
                                   seg_file='sdssj0851+3331-f160w_seg.fits',
                                   catalog='sdssj0851+3331-f160w.cat',
                                   cpu_count=4)
    
    beams = grp.get_beams(LENS_SEG_ID, size=80, beam_id='A')
    
    mb = grizli.multifit.MultiBeam(beams, fcontam=0., group_name='sgas-lens')
    
    # Save them to a file for quick I/O
    mb.write_beam_fits()
    
    # id661acmq_flt_00253.g141.A.fits
    # id661acnq_flt_00253.g141.A.fits
    # id661bcuq_flt_00253.g141.A.fits
    # id661bcvq_flt_00253.g141.A.fits
    # etc.

def load_beams():
    """
    Load the beams from the saved files
    """
    import glob
    import grizli
    
    beam_files = glob.glob('*{0}.g141.A.fits'.format(LENS_SEG_ID))
    # Load from FITS files
    mb = grizli.multifit.MultiBeam(beam_files, fcontam=0., group_name='sgas-lens')
    
    return mb
    
def make_component_segments(mb, region_file='sgas-lens-components.reg', force_recompute=False):
    """
    Make custom segmentation file for each beam in the `mb` list.
    
    The `region_file` file is in the sgas-lens/sgas/data repository 
    """
    import os
    import numpy as np
    
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import pyregion
    
    import sgas
    
    # Load the segment file if already exists and if not forced to recompute
    if os.path.exists('component_segments.npy') & (not force_recompute):
        component_segments, old_region_file = np.load('component_segments.npy')
        if (old_region_file == region_file):
            for ib, beam in enumerate(mb.beams):
                beam.beam.seg = component_segments[ib]
            
            return mb
    
    # Read the regions file    
    regions = pyregion.open(os.path.join(sgas.get_data_path(), region_file))
    Nreg = len(regions)
    
    component_segments = []
    
    # Loop through beams, make segments for each one of them
    for beam in mb.beams:
        regs = regions.as_imagecoord(beam.direct.header)
        
        #initialize a new segmentation array
        seg = beam.beam.seg*0
        
        # update with region segments
        for i in range(Nreg):
            mask_i = regs[i:i+1].get_mask(header=beam.direct.header, shape=beam.beam.sh)
            
            # update nonzero segment pixels
            clip = seg == 0
            seg[clip] += mask_i[clip]*(i+1)
        
        beam.beam.seg = seg
        component_segments.append(seg)
    
    np.save('component_segments.npy', [component_segments, region_file])
    return mb
    
def test():
    # SHow an example
    model = mb.beams[0].beam.model*0.
    for i in range(6):
        set_id(mb, id=i+1)
        mb.compute_model()
        model += mb.beams[0].beam.model
        ds9.view(model)
        
def set_id(mb, id=1):
    # Set the id on the MultiBeam object for a given component
    mb.id = id
    for beam in mb.beams:
        beam.id = id
        beam.beam.id = id
        
        ### Initialize for fits
        beam.flat_flam = beam.compute_model(in_place=False)
        
        ### OK data where the 2D model has non-zero flux
        beam.fit_mask = (~beam.mask.flatten()) & (beam.ivar.flatten() != 0)
        beam.fit_mask &= (beam.flat_flam > 0.01*beam.flat_flam.max())
        
        try:
            delattr(beam, 'optimal_profile')
        except:
            pass
    
    mb.flat_flam = np.hstack([b.flat_flam for b in mb.beams])                            
    mb.fit_mask = np.hstack([b.fit_mask for b in mb.beams])
    mb.DoF = mb.fit_mask.sum()
    
    