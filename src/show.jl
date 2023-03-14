#Common text for display multivectors.

Base.show(io::IO, ::MIME"text/plain", mv::Even) = print(io, "", mvtype(mv))

function Base.show(io::IO, ::MIME"text/plain", mvs::Vector{Even})
    n= length(mvs)
    println(io,n,"-element Vector{MVeven}")
    for i in eachindex(mvs)
    println(io, " ", mvtype(mvs[i]))
    end
end


Base.show(io::IO, ::MIME"text/plain", mv::Odd) = print(io, "", mvtype(mv))

function Base.show(io::IO, ::MIME"text/plain", mvs::Vector{Odd})
    n= length(mvs)
    println(io,n,"-element Vector{MVodd}")
    for i in eachindex(mvs)
    println(io, " ", mvtype(mvs[i]))
    end
end

