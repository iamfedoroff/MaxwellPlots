module MaxwellPlots

import FFTW
import GLMakie as mak
import HDF5
import KernelAbstractions: @index, @kernel, get_backend
import Printf: @sprintf

export plot_geometry,
       plot_waveform,
       plot1D,
       inspect1D,
       plot1D_line,
       plot2D,
       inspect2D,
       inspect2D_xsec,
       inspect2D_poynting,
       plot2D_poynting_averaged,
       inspect3D,
       inspect3D_xsec,
       inspect3D_volume,
       plot3D,
       plot3D_diff,
       plot3D_xsec,
       plot3D_xsec_diff,
       plot_viewpoints,
       plot_viewpoints_spectrum,
       plot_viewpoints_polarization,
       rough_plot

include("util.jl")
include("recipes3D.jl")
include("plot_geometry.jl")
include("plots1D.jl")
include("plots2D.jl")
include("plots3D.jl")
include("viewpoints.jl")
include("rough.jl")

end
