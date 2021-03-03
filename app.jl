using Dash, DashHtmlComponents, DashCoreComponents
using DataFrames, UrlDownload 
using LiteQTL

############################ Structures #################################
mutable struct ScanArgs 
    genofile::String
    phenofile::String
    covarfile::String 
    exportMatrix::Bool
    usegpu::Bool
    ScanArgs(genofile, phenofile, covarfile, exportMatrix, usegpu) = new(genofile, phenofile,covarfile,exportMatrix, usegpu)
end

############################# functions ####################################

function getUserInput(genofile, phenofile, covarfile, useroptions)

    args = ScanArgs("", "", "", false, false)

    args.genofile = genofile
    args.phenofile = phenofile
    args.covarfile = covarfile

    if isnothing(useroptions)  
        println("is nothing or 0  ")
    else
        println("else $useroptions")
        for i in 1:length(useroptions)
            if useroptions[i] == "export_matrix"
                args.exportMatrix = true
            elseif useroptions[i] == "usegpu"
                args.usegpu = true 
            end
        end
    end
    return args
end

function callScan(args::ScanArgs)
    G = get_geno_data(args.genofile, Float64)
    Y = get_pheno_data(args.phenofile, Float64, transposed=false)
    n = size(Y,1)
    m = size(Y,2)
    p = size(G,2)
    return LiteQTL.scan(Y, G, n; export_matrix = args.exportMatrix, usegpu = args.usegpu);
end

############################# Initializations ##############################
app = dash()

dropdown_options = [
    Dict("label" => "Export Matrix", "value" => "export_matrix"),
    Dict("label" => "Use GPU", "value" => "usegpu"),
]

############################# Run ##########################################
app.layout = html_div() do
    html_h1("LiteQTL.jl Demonstration"),

    html_h4("LiteQTL: A lightweight Julia package for eQTL genome scans near real-time."),

    dcc_input(id = "inputgenofile", value = "Genotype file location", type = "text"),

    dcc_input(id = "inputphenofile", value = "Phenotype file location", type = "text"),

    dcc_input(id = "inputcovarfile", value = "Covariates file location", type = "text"),
    
    dcc_checklist(id = "inputoptions", options = dropdown_options), #, value = ["usegpu"]

    html_button(id = "run-button", children = "run", n_clicks = 0),
    
    html_div(id = "outputlod")
    
end

callback!(app, 
    Output("outputlod", "children"),
    Input("run-button", "n_clicks"),
    State("inputgenofile", "value"),
    State("inputphenofile", "value"),
    State("inputcovarfile", "value"),
    State("inputoptions", "value")
    
) do clicks, genofile, phenofile, covarfile, useroptions
    if clicks >= 1
        inputargs = getUserInput(genofile, phenofile, covarfile, useroptions)
        lod = callScan(inputargs)
    
        return display(lod[10, :])
    end
end

run_server(app, "0.0.0.0", 5240, debug = true)

# ../LiteQTL.jl/data/processed/spleen-bxd-genoprob.csv
# ../LiteQTL.jl/data/processed/spleen-pheno-nomissing.csv
