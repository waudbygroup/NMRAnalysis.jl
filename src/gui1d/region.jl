function Region(x1, x2, label::String, expt::Experiment)
    x1 = Observable(min(x1,x2))
    x2 = Observable(max(x1,x2))
    label = Observable(label)
    integrals = @lift integrate($x1, $x2, expt)
    fitparameters = @lift fit($integrals, expt)
    
    Region(x1, x2, label, integrals, fitparameters)
end