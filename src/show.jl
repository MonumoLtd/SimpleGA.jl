# Common implementation for displaying multivectors.
# This is intended to be included in the various gaXXX.jl files, where the
#   `Even` and `Odd` names will refer to types in that particular module.
Base.show(io::IO, mv::Even) = print(io, mv_to_text(mv))
Base.show(io::IO, mv::Odd) = print(io, mv_to_text(mv))
