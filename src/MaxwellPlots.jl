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
       plot2D,
       inspect2D,
       inspect2D_xsec,
       inspect2D_poynting,
       plot2D_poynting_averaged,
       inspect3D,
       inspect3D_xsec,
       inspect3D_poynting,
       inspect3D_poynting_xsec,
       inspect3D_poynting_averaged,
       plot3D_poynting_averaged,
       plot3D_poynting_averaged_diff,
       plot3D_poynting_averaged_xsec,
       plot3D_poynting_averaged_xsec_diff,
       plot_viewpoints,
       plot_viewpoints_spectrum,
       plot_viewpoints_polarization

include("util.jl")
include("recipes3D.jl")
include("plots1D.jl")
include("plots2D.jl")
include("plots3D.jl")
include("viewpoints.jl")

end
