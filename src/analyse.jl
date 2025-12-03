function analyse(filename)
    # generic analysis based on annotations and experiment types
    types, features = types_and_features(filename)
    if isempty(types)
        @warn "Could not determine experiment type for file: $filename. Unable to proceed with analysis."
        return
    end
    if "1d" in types
        if "calibration" in types
            return analyse_1d_calibration(filename, types, features)
        end
    end
    @info "No automatic analysis routine available for experiment types: $(join(types, ", "))"
end

function analyse_1d_calibration(filename, types, features)
    if "1d" in types && "nutation" in features
        @info "Analysing calibration by nutation on $filename"
        analyse_1d_nutation(filename)
    end
end
