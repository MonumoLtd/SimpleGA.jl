#Common text for display multivectors.

Base.show(io::IO, mv::Even) = print(io, mvtype(mv))
Base.show(io::IO, mv::Odd) = print(io, mvtype(mv))
