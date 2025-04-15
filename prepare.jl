sing Base.Threads

include("batch/1cusp_manifest.txt")

@threads for l in isosigs
    println(l)
    run(`python3 prepare.py $(l)`)
end
