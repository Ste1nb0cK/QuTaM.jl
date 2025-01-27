using BackAction
using Documenter
using DocumenterCitations
# Set up the bibliography
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

# GPT's explanation:
# This sets up metadata for BackAction with the :DocTestSetup key.
# The :DocTestSetup specifies that when running documentation tests (doctests),
# the statement using BackAction should be included.
# This ensures the module is loaded before running any test snippets in the documentation.
# recursive = true applies this metadata recursively to submodules, ensuring all of
# BackAction's components are included.
# DocMeta.setdocmeta!(BackAction, :DocTestSetup, :(using BackAction); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers


makedocs(;
  modules = [BackAction], # Specify the module to document
  doctest = true, # Run tests in the documentation to ensure the examples work
  linkcheck = false, # Disable link checking
  authors = "Nicolás A. Niño <ninino@unal.edu.co>",
  repo = "https://github.com/Ste1nb0cK/BackAction.jl/blob/{commit}{path}#{line}",
  sitename = "BackAction.jl",
  # Configure output format
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") ==  true,
    # canonical = "https://abelsiqueira.github.io/COPIERTemplate.jl", # Canonical URL

  ),
  # Define the structure of the site
  pages = [
    "Home" => "index.md",
      "Getting Started" => "getting_started.md",
      "System and SimulParameters" => "system_simulparams.md",
      "Jump Unraveling" => "jump_unraveling.md",
      "Monitoring Quantum Metrology" => "monitoring_metrology.md",
      "Bibliography" => "bibliography.md"
  ],
  plugins=[bib]
)

deploydocs(; repo = "github.com/Ste1nb0cK/BackAction.jl"
           )
