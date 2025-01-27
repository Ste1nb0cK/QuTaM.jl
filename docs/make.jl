using BackAction
using Documenter
using DocumenterCitations
# Set up the bibliography
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))


makedocs(;
  modules = [BackAction], # Specify the module to document
  authors = "Nicolás A. Niño <ninino@unal.edu.co>",
  sitename = "BackAction.jl",
  # Configure output format
  format = Documenter.HTML(repolink = "..."),

  # Define the structure of the site
  pages = [
    "Home" => "index.md",
      "Getting Started" => "getting_started.md",
      "System and SimulParameters" => "system_simulparams.md",
      "Jump Unraveling" => "jump_unraveling.md",
      "Monitoring Quantum Metrology" => "monitoring_metrology.md",
      "Bibliography" => "bibliography.md"
  ],
  repo = "https://github.com/Ste1nb0cK/BackAction.jl/blob/{commit}{path}#{line}",
  plugins=[bib]
)

deploydocs(;
           repo = "github.com/Ste1nb0cK/BackAction.jl",
           devurl="dev")
