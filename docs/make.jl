using BackAction
using Documenter
using DocumenterCitations
# Set up the bibliography
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
  modules = [BackAction], # Specify the module to document
  doctest = true,
  authors = "Nicolás A. Niño <ninino@unal.edu.co>",
  sitename = "BackAction.jl",
  repo = "https://github.com/Ste1nb0cK/BackAction.jl/blob/{commit}{path}#{line}",
  # Configure output format
  format = Documenter.HTML(;
                           prettyurls = true,
                           canonical = "https://github.com/Ste1nb0cK/BackAction.jl",
                           assets = ["assets/style.css"],
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

deploydocs(;
           repo = "github.com/Ste1nb0cK/BackAction.jl",
           push_preview=true)
