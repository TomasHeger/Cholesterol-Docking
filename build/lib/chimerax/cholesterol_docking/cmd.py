from chimerax.core.commands import CmdDesc, ModelsArg, StringArg, FloatArg, IntArg
from .engine import run_docking

def choldock7(session, map_models, *, protein_models, out='cholesterol_v7', detection_level=.141,
              support_level=.05, orientations=1500, maxsites=60, max_tilt=45., mask_radius=4.,
              protein_map_resolution=6., blob_smoothing=3., blob_threshold_fraction=.10,
              protein_mask_radius=2.5, seed=12345):
    return run_docking(session,map_models,protein_models,out,
        detection_level=detection_level,support_level=support_level,orientations=orientations,
        maxsites=maxsites,max_tilt=max_tilt,mask_radius=mask_radius,
        protein_map_resolution=protein_map_resolution,blob_smoothing=blob_smoothing,
        blob_threshold_fraction=blob_threshold_fraction,protein_mask_radius=protein_mask_radius,seed=seed)

choldock7_desc=CmdDesc(
    required=[('map_models',ModelsArg)],
    keyword=[('protein_models',ModelsArg),('out',StringArg),('detection_level',FloatArg),
             ('support_level',FloatArg),('orientations',IntArg),('maxsites',IntArg),
             ('max_tilt',FloatArg),('mask_radius',FloatArg),('protein_map_resolution',FloatArg),
             ('blob_smoothing',FloatArg),('blob_threshold_fraction',FloatArg),
             ('protein_mask_radius',FloatArg),('seed',IntArg)],
    required_arguments=['protein_models'],
    synopsis='Automated multi-site cholesterol docking into cryo-EM density')
