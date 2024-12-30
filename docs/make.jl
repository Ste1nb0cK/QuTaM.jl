using QuTaM
using Documenter

DocMeta.setdocmeta!(QuTaM, :DocTestSetup, :(using QuTaM); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

function nice_name(file)
  file = replace(file, r"^[0-9]*-" => "")
  if haskey(page_rename, file)
    return page_rename[file]
  end
  return splitext(file)[1] |> x -> replace(x, "-" => " ") |> titlecase
end

makedocs(;
  modules = [QuTaM],
  doctest = true,
  linkcheck = false, # Rely on Lint.yml/lychee for the links
  authors = "Nicolás A. Niño <ninino@unal.edu.co>",
  repo = "https://github.com/Ste1nb0cK/QuTaM.jl/blob/{commit}{path}#{line}",
  sitename = "QuTaM.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") ==  true,
    canonical = "https://abelsiqueira.github.io/COPIERTemplate.jl"
  ),
  pages = [
    "Home" => "index.md", "Tutorial" => "tutorial.md",
      "System and SimulParameters" => "system_simulparams.md", "Reference" => "reference.md",
    # [
    #   nice_name(file) => file for
    #   file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && splitext(file)[2] == ".md"
    # ]
  ],
)

deploydocs(; repo = "github.com/Ste1nb0cK/QuTaM.jl", push_preview = false)
