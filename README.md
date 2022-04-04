# Template for JuDO Packages

## Create GitHub Repository
1. Start by clicking the green **'Use this template'** button;
2. Name your repository with the `.jl` suffix, like `Pizza.jl`;
3. Decide on its visibility: `Public`/`Private`;
4. Leave **'Include all branches'** unticked;
5. Clone the repository to your machine* (e.g., inside `\Documents`).

## Generate Julia Package
1. Open Julia and create a package template by running:
```julia 
julia> using Pkg; Pkg.add("PkgTemplates");
julia> using PkgTemplates
julia> t = Template(
    user="JuDO-dev",
    dir=pwd(),
    julia=v"1.6",
    plugins=[
        !License,
        Git(branch="dev"),
        GitHubActions(linux=true, x64=true, x86=true, extra_versions=[v"1.7", "nightly"]),
        Codecov(),
        Documenter{GitHubActions}()])
```
2. Generate the package files by running:
```julia
julia> t("Pizza") # NB: without the ".jl" suffix
```

## Assemble Julia Package
1. Copy the contents of `Pizza` into `Pizza.jl`, overwritting `README.md`;
2. Commit changes and push to origin*.

*You may use [GitHub Desktop](https://desktop.github.com/).
