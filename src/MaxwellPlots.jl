module MaxwellPlots

import FFTW
import GLMakie as mak
import HDF5
import KernelAbstractions: @index, @kernel, get_backend
import Printf: @sprintf

export plot_geometry,
       plot_waveform,
       inspect,
       inspect1D,
       inspect2D,
       inspect3D,
       plot1D,
       plot1D_line,
       plot2D,
       inspect2D_xsec,
       inspect2D_poynting,
       plot2D_poynting_averaged,
       inspect3D_xsec,
       inspect3D_volume,
       plot3D,
       plot3D_diff,
       plot3D_xsec,
       plot3D_xsec_diff,
       plot_monitors,
       plot_monitors_spectrum,
       plot_monitors_polarization,
       plot_monitors_polarization_xsec,
       rough_plot,
       plot_spatial_spectrum,
       inspect_spatial_spectrum,
       plot_integrated_spatial_spectrum

include("util.jl")
include("inspect.jl")
include("recipes3D.jl")
include("plot_geometry.jl")
include("plots1D.jl")
include("plots2D.jl")
include("plots3D.jl")
include("monitors.jl")
include("rough.jl")
include("spatial_spectra.jl")

end
