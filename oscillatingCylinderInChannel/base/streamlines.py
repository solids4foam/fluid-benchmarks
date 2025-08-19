# para-anim.py  —  headless-friendly velocity+black-streamlines animation
# Run with:
#   pvbatch --opengl-window-backend=EGL para-anim.py

from paraview.simple import *
import os

# --------- inputs ---------
case_path      = "./case.foam"
out_movie      = "animation.avi"   # change to .avi if your build lacks mp4
fps            = 24
seed_count     = 50               # streamline density
line_width     = 2.0
start_time_index = 0
# --------------------------

paraview.simple._DisableFirstRenderCameraReset()

# Reader (2D front slice as in your script)
print(f"Reading {case_path}")
src = OpenFOAMReader(FileName=case_path)
src.Createcelltopointfiltereddata = True

# Keep only the front patch if your case is extruded in Z:
src.MeshRegions = ['patch/front']
src.UpdatePipeline()

# Animation scene follows data timesteps
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
timesteps = list(scene.TimeKeeper.TimestepValues or [0.0])

# View
view = GetActiveViewOrCreate('RenderView')
view.OrientationAxesVisibility = 1
view.UseColorPaletteForBackground = 0
view.Background = [1, 1, 1]                 # white background
view.InteractionMode = '2D'
view.CameraParallelProjection = 1
view.OrientationAxesVisibility = 0   # hide triad
view.CenterAxesVisibility = 0
view.ResetCamera(False, 0.99)
view.CameraParallelScale *= 0.70     # ~30% zoom-in; tweak 0.5–0.8 to get the zoom right

# --- base coloured field (U magnitude) ---
disp = Show(src, view, 'GeometryRepresentation')

# Colour by velocity or just white
ColorBy(disp, ('POINTS', 'U', 'Magnitude'))
#ColorBy(disp, None)              # disable scalar colouring
#disp.DiffuseColor = [1.0, 1.0, 1.0]   # RGB white

disp.RescaleTransferFunctionToDataRange(False, True)  # rescale over ALL timesteps

# Display mesh edges
#edges = ExtractEdges(Input=src)
#edgesD = Show(edges, view)
#edgesD.DiffuseColor = [0, 0, 0]
#edgesD.LineWidth = 1.0
#ColorBy(edgesD, None)   # no scalar colouring

uLUT = GetColorTransferFunction('U')
uPWF = GetOpacityTransferFunction('U')
uLUT = GetColorTransferFunction('U')
uLUT.Discretize = 1
uLUT.NumberOfTableValues = 12
# Optional: pick a preset, then keep it discrete
try:
    #uLUT.ApplyPreset('Cool to Warm (Extended)', True)
    uLUT.ApplyPreset('Warm to Cool (Extended)', True)
    #uLUT.ApplyPreset('Viridis (matplotlib)', True)
    #uLUT.InvertTransferFunction = 1 
except:
    pass

view.OrientationAxesVisibility = 0   # hide bottom-left triad
view.AxesGrid.Visibility = 0         # grid off (usually default)
# If a 2D axis annotation is on:
try:
    view.UseSeparateRenderers = 0
except:
    pass

# --- black outline of the geometry
edges = FeatureEdges(Input=src)
edgesD = Show(edges, view, 'GeometryRepresentation')
ColorBy(edgesD, None)
edgesD.DiffuseColor = [0, 0, 0]
edgesD.LineWidth = 10.0


print("Creating streamlines")
# create a new 'Evenly Spaced Streamlines 2D'
stream = EvenlySpacedStreamlines2D(registrationName='EvenlySpacedStreamlines2D1', Input=src)
stream.SeparatingDistance = 1.0 #1.5
#stream.SeparatingDistanceRatio = 1
stream.StartPosition = [0.5, 0.0, 0.0]
stream.MaximumSteps = 100000
stream.ComputeVorticity = False
stream.ClosedLoopMaximumDistance = 0.1

streamD = Show(stream, view, 'GeometryRepresentation')
ColorBy(streamD, None)
streamD.DiffuseColor = [1, 1, 1] # white
#streamD.DiffuseColor = [0, 0, 0] # black
streamD.Opacity = 0.2
streamD.LineWidth = line_width
streamD.RenderLinesAsTubes = 0  # set 1 for thicker tube look (slower)


# --- Tight 2D framing (no big white borders) ---
view = GetActiveViewOrCreate('RenderView')
view.InteractionMode = '2D'
view.CameraParallelProjection = 1
view.OrientationAxesVisibility = 0
view.CenterAxesVisibility = 0
view.UseColorPaletteForBackground = 0
view.Background = [1,1,1]

# Bounds and aspect
xmin,xmax,ymin,ymax,_,_ = GetActiveSource().GetDataInformation().GetBounds()
cx = 0.5*(xmin+xmax); cy = 0.5*(ymin+ymax)
Wworld = (xmax - xmin); Hworld = (ymax - ymin)
aspect = Wworld / Hworld if Hworld > 0 else 1.0

# Camera exactly on the slice, tiny padding
view.CameraFocalPoint = [cx, cy, 0.0]
view.CameraPosition   = [cx, cy, 1.0]        # any z>0 is fine in 2D
view.CameraViewUp     = [0, 1, 0]
view.CameraParallelScale = 0.5*Hworld * 1.02 # ~2% border

# Choose output resolution to match domain aspect ratio (prevents letterboxing)
Hpx = 1080                                   # pick your height
Wpx = int(round(Hpx * aspect))
image_res = [Wpx, Hpx]
print("ImageResolution:", image_res)


print(f"Exporting {out_movie}")
try:
    SaveAnimation(out_movie, view,
                  ImageResolution=image_res,
                  FrameRate=fps,
                  FrameWindow=[start_time_index, len(timesteps)-1])
    print(f"Wrote {out_movie}")
except Exception as e:
    # Fallback to PNG frames if codec not available
    print(f"SaveAnimation to movie failed ({e}); writing PNG frames instead…")
    SaveAnimation("frame.png", view,
                  ImageResolution=image_res,
                  FrameRate=fps,
                  FrameWindow=[start_time_index, len(timesteps)-1])
    print("Wrote PNG frames: frame.0000.png …")

