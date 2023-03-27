#Common text for display multivectors.

Base.show(io::IO, mv::Even) = print(io, mv_to_text(mv))
Base.show(io::IO, mv::Odd) = print(io, mv_to_text(mv))
