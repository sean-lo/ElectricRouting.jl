using Pkg
Pkg.activate(".")
include("utils.jl")

abstract type SubpathLabel <: Label end

function get_subpath_vkey_type(
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
)
    """
    (
        cost::Float64,
        if (
            customer_service isa NgRoute
            && cuts isa SR3Cuts
        )
            cut_flabels::BitVector (n_cuts,)
        elseif (
            customer_service isa NgRoute
            && cuts isa LmSR3Cuts
        )
            cut_flabels::BitVector (n_cuts,)
            cut_blabels::BitVector (n_cuts,)
            cut_mlabels::BitVector (n_cuts,)
        end
        if use_time_windows isa YesTimeWindows
            time_end_earliest::Int,
            time_start_latest::Int,
            max(time_end_earliest - time_start_latest, time_taken)::Int,
        else
            time_taken::Int,
        end
        charge::Int,
        if use_load isa YesLoad
            load::Int,
        end
        if customer_service isa Elementary
            served::BitVector (n_customers,)
        elseif customer_service isa NgRoute
            ng_fset::BitVector (n_nodes,)
            ng_residue::BitVector (n_nodes,)
            ng_bset::BitVector (n_nodes,)
        end
    )
    """
    if use_load isa YesLoad
        if use_time_windows isa YesTimeWindows
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, Int, Int, 
                    Int, Int, 
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, Int, Int, 
                    Int, Int, 
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64,     
                        Int, Int, Int, 
                        Int, Int, 
                        BitVector, BitVector, BitVector,
                    }
                elseif cuts isa SR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector,
                        Int, Int, Int, 
                        Int, Int, 
                        BitVector, BitVector, BitVector, 
                    }
                elseif cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector, BitVector, BitVector,
                        Int, Int, Int, 
                        Int, Int, 
                        BitVector, BitVector, BitVector, 
                    }
                end
            end
        elseif use_time_windows isa NoTimeWindows
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, 
                    Int, Int,
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, 
                    Int, Int,
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, 
                        Int, Int, 
                        BitVector, BitVector, BitVector,
                    }
                elseif cuts isa SR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector,
                        Int, 
                        Int, Int, 
                        BitVector, BitVector, BitVector, 
                    }
                elseif cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector, BitVector, BitVector,
                        Int, 
                        Int, Int, 
                        BitVector, BitVector, BitVector, 
                    }
                end
            end
        end
    elseif use_load isa NoLoad
        if use_time_windows isa YesTimeWindows
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, Int, Int,
                    Int,
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, Int, Int,
                    Int,
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, Int, Int,
                        Int,
                        BitVector, BitVector, BitVector,
                    }
                elseif cuts isa SR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector,
                        Int, Int, Int,
                        Int,
                        BitVector, BitVector, BitVector, 
                    }
                elseif cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector, BitVector, BitVector,
                        Int, Int, Int,
                        Int,
                        BitVector, BitVector, BitVector, 
                    }
                end
            end
        elseif use_time_windows isa NoTimeWindows
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, 
                    Int,
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, 
                    Int,
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, 
                        Int, 
                        BitVector, BitVector, BitVector,
                    }
                elseif cuts isa SR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector,
                        Int, 
                        Int, 
                        BitVector, BitVector, BitVector, 
                    }
                elseif cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector, BitVector, BitVector,
                        Int, 
                        Int, 
                        BitVector, BitVector, BitVector, 
                    }
                end
            end
        end
    end
end

function get_subpath_vkey_fkey_type(
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
)
    return Tuple{
        get_subpath_vkey_type(
            use_load, 
            use_time_windows,
            customer_service,
            cuts,
        ).parameters...,
        # starting_node, current_node,
        Int, Int,
    }
end

function get_path_vkey_type(
    use_load::Load,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
)
    """
    (
        cost::Float64,
        if (
            customer_service isa NgRoute
            && (cuts isa SR3Cuts || cuts isa LmSR3Cuts)
        )
            cut_flabels::BitVector (n_cuts,)
        end
        time::Int,
        charge::Int,
        if use_load isa YesLoad
            load::Int,
        end
        if charge_costs isa HetCharge
            charge_rebalance_slacks::Vector{Int} (charge_cost_nlevels - 1,)
        end
        if customer_service isa Elementary
            served::BitVector (n_customers,)
        elseif customer_service isa NgRoute
            ng_fset::BitVector (n_nodes,)
        end
    )
    """
    if use_load isa YesLoad
        if charge_costs isa HomCharge
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, Int, Int,
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, Int, Int,
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, Int, Int, 
                        BitVector, 
                    }
                elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector,
                        Int, Int, Int, 
                        BitVector, 
                    }
                end
            end
        elseif charge_costs isa HetCharge
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, Int, Int,
                    Vector{Int},
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, Int, Int,
                    Vector{Int},
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, Int, Int, 
                        Vector{Int},
                        BitVector, 
                    }
                elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector, 
                        Int, Int, Int, 
                        Vector{Int},
                        BitVector, 
                    }
                end
            end
        end
    elseif use_load isa NoLoad
        if charge_costs isa HomCharge
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, Int,
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, Int,
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, Int, 
                        BitVector, 
                    }
                elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector,
                        Int, Int, 
                        BitVector, 
                    }
                end
            end
        elseif charge_costs isa HetCharge
            if customer_service isa Elementary
                return Tuple{
                    Float64, 
                    Int, Int,
                    Vector{Int},
                    BitVector,
                }
            elseif customer_service isa NoService
                return Tuple{
                    Float64, 
                    Int, Int,
                    Vector{Int},
                }
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return Tuple{
                        Float64, 
                        Int, Int, 
                        Vector{Int},
                        BitVector, 
                    }
                elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                    return Tuple{
                        Float64, 
                        BitVector, 
                        Int, Int, 
                        Vector{Int},
                        BitVector, 
                    }
                end
            end
        end
    end
end

function get_path_vkey_fkey_type(
    use_load::Load,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
)
    return Tuple{
        get_path_vkey_type(
            use_load, 
            charge_costs,
            customer_service,
            cuts,
        ).parameters...,
        # starting_node, current_node,
        Int, Int,
    }
end

abstract type YesLoadSubpathLabel <: SubpathLabel end
abstract type YesLoadYesTimeWindowsSubpathLabel <: YesLoadSubpathLabel end

mutable struct YesLoadYesTimeWindowsElementarySubpathLabel <: SubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    served::BitVector
end
get_vkey(s::YesLoadYesTimeWindowsElementarySubpathLabel) = (
    s.cost,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.load,
    s.served,
)
mutable struct YesLoadYesTimeWindowsNoServiceSubpathLabel <: SubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
end
get_vkey(s::YesLoadYesTimeWindowsNoServiceSubpathLabel) = (
    s.cost,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.load,
)

abstract type YesLoadYesTimeWindowsNgRouteSubpathLabel <: YesLoadYesTimeWindowsSubpathLabel end

mutable struct YesLoadYesTimeWindowsNgRouteNoCutsSubpathLabel <: YesLoadYesTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
end
get_vkey(s::YesLoadYesTimeWindowsNgRouteNoCutsSubpathLabel) = (
    s.cost,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.load,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct YesLoadYesTimeWindowsNgRouteSR3CutsSubpathLabel <: YesLoadYesTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
end
get_vkey(s::YesLoadYesTimeWindowsNgRouteSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.load,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct YesLoadYesTimeWindowsNgRouteLmSR3CutsSubpathLabel <: YesLoadYesTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
    cut_blabels::BitVector
    cut_mlabels::BitVector
end
get_vkey(s::YesLoadYesTimeWindowsNgRouteLmSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.cut_blabels,
    s.cut_mlabels,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.load,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)


abstract type YesLoadNoTimeWindowsSubpathLabel <: YesLoadSubpathLabel end
mutable struct YesLoadNoTimeWindowsElementarySubpathLabel <: YesLoadNoTimeWindowsSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    served::BitVector
end
get_vkey(s::YesLoadNoTimeWindowsElementarySubpathLabel) = (
    s.cost,
    s.time_taken,
    s.charge_taken,
    s.load,
    s.served,
)
mutable struct YesLoadNoTimeWindowsNoServiceSubpathLabel <: YesLoadNoTimeWindowsSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
end
get_vkey(s::YesLoadNoTimeWindowsNoServiceSubpathLabel) = (
    s.cost,
    s.time_taken,
    s.charge_taken,
    s.load,
)

abstract type YesLoadNoTimeWindowsNgRouteSubpathLabel <: YesLoadNoTimeWindowsSubpathLabel end

mutable struct YesLoadNoTimeWindowsNgRouteNoCutsSubpathLabel <: YesLoadNoTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
end
get_vkey(s::YesLoadNoTimeWindowsNgRouteNoCutsSubpathLabel) = (
    s.cost,
    s.time_taken,
    s.charge_taken,
    s.load,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct YesLoadNoTimeWindowsNgRouteSR3CutsSubpathLabel <: YesLoadNoTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
end
get_vkey(s::YesLoadNoTimeWindowsNgRouteSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.time_taken,
    s.charge_taken,
    s.load,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct YesLoadNoTimeWindowsNgRouteLmSR3CutsSubpathLabel <: YesLoadNoTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    load::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
    cut_blabels::BitVector
    cut_mlabels::BitVector
end
get_vkey(s::YesLoadNoTimeWindowsNgRouteLmSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.cut_blabels,
    s.cut_mlabels,
    s.time_taken,
    s.charge_taken,
    s.load,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

abstract type NoLoadSubpathLabel <: SubpathLabel end

abstract type NoLoadYesTimeWindowsSubpathLabel <: NoLoadSubpathLabel end
mutable struct NoLoadYesTimeWindowsElementarySubpathLabel <: NoLoadYesTimeWindowsSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    nodes::Vector{Int}
    served::BitVector
end
get_vkey(s::NoLoadYesTimeWindowsElementarySubpathLabel) = (
    s.cost,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.served,
)
mutable struct NoLoadYesTimeWindowsNoServiceSubpathLabel <: NoLoadYesTimeWindowsSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    nodes::Vector{Int}
end
get_vkey(s::NoLoadYesTimeWindowsNoServiceSubpathLabel) = (
    s.cost,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
)

abstract type NoLoadYesTimeWindowsNgRouteSubpathLabel <: NoLoadYesTimeWindowsSubpathLabel end

mutable struct NoLoadYesTimeWindowsNgRouteNoCutsSubpathLabel <: NoLoadYesTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
end
get_vkey(s::NoLoadYesTimeWindowsNgRouteNoCutsSubpathLabel) = (
    s.cost,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct NoLoadYesTimeWindowsNgRouteSR3CutsSubpathLabel <: NoLoadYesTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
end
get_vkey(s::NoLoadYesTimeWindowsNgRouteSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct NoLoadYesTimeWindowsNgRouteLmSR3CutsSubpathLabel <: NoLoadYesTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    time_end_earliest::Int
    time_start_latest::Int
    charge_taken::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
    cut_blabels::BitVector
    cut_mlabels::BitVector
end
get_vkey(s::NoLoadYesTimeWindowsNgRouteLmSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.cut_blabels,
    s.cut_mlabels,
    s.time_end_earliest,
    - s.time_start_latest,
    max(s.time_end_earliest - s.time_start_latest, s.time_taken),
    s.charge_taken,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)


abstract type NoLoadNoTimeWindowsSubpathLabel <: NoLoadSubpathLabel end

mutable struct NoLoadNoTimeWindowsElementarySubpathLabel <: NoLoadNoTimeWindowsSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    nodes::Vector{Int}
    served::BitVector
end
get_vkey(s::NoLoadNoTimeWindowsElementarySubpathLabel) = (
    s.cost,
    s.time_taken,
    s.charge_taken,
    s.served,
)
mutable struct NoLoadNoTimeWindowsNoServiceSubpathLabel <: NoLoadNoTimeWindowsSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    nodes::Vector{Int}
end
get_vkey(s::NoLoadNoTimeWindowsNoServiceSubpathLabel) = (
    s.cost,
    s.time_taken,
    s.charge_taken,
)

abstract type NoLoadNoTimeWindowsNgRouteSubpathLabel <: NoLoadNoTimeWindowsSubpathLabel end

mutable struct NoLoadNoTimeWindowsNgRouteNoCutsSubpathLabel <: NoLoadNoTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
end
get_vkey(s::NoLoadNoTimeWindowsNgRouteNoCutsSubpathLabel) = (
    s.cost,
    s.time_taken,
    s.charge_taken,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct NoLoadNoTimeWindowsNgRouteSR3CutsSubpathLabel <: NoLoadNoTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
end
get_vkey(s::NoLoadNoTimeWindowsNgRouteSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.time_taken,
    s.charge_taken,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)

mutable struct NoLoadNoTimeWindowsNgRouteLmSR3CutsSubpathLabel <: NoLoadNoTimeWindowsNgRouteSubpathLabel
    cost::Float64
    time_taken::Int
    charge_taken::Int
    nodes::Vector{Int}
    ng_fset::BitVector
    ng_residue::BitVector
    ng_bset::BitVector
    cut_flabels::BitVector
    cut_blabels::BitVector
    cut_mlabels::BitVector
end
get_vkey(s::NoLoadNoTimeWindowsNgRouteLmSR3CutsSubpathLabel) = (
    s.cost,
    s.cut_flabels,
    s.cut_blabels,
    s.cut_mlabels,
    s.time_taken,
    s.charge_taken,
    s.ng_fset,
    s.ng_residue,
    s.ng_bset,
)


function get_subpath_label_type(
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
)
    if use_load isa YesLoad
        if customer_service isa Elementary
            return YesLoadElementarySubpathLabel
        elseif customer_service isa NoService
            return YesLoadNoServiceSubpathLabel
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return YesLoadNgRouteNoCutsSubpathLabel
            elseif cuts isa SR3Cuts
                return YesLoadNgRouteSR3CutsSubpathLabel
            elseif cuts isa LmSR3Cuts
                return YesLoadNgRouteLmSR3CutsSubpathLabel
            end
        end
    elseif use_load isa NoLoad
        if customer_service isa Elementary
            return NoLoadElementarySubpathLabel
        elseif customer_service isa NoService
            return NoLoadNoServiceSubpathLabel
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return NoLoadNgRouteNoCutsSubpathLabel
            elseif cuts isa SR3Cuts
                return NoLoadNgRouteSR3CutsSubpathLabel
            elseif cuts isa LmSR3Cuts
                return NoLoadNgRouteLmSR3CutsSubpathLabel
            end
        end
    end
end

abstract type PPathLabel <: Label end
abstract type YesLoadPPathLabel <: PPathLabel end
abstract type YesLoadHomChargePPathLabel <: YesLoadPPathLabel end

mutable struct YesLoadHomChargeElementaryPPathLabel <: YesLoadHomChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadElementarySubpathLabel} # TODO: replace with dict keys to nondominated subpaths
    charging_actions::Vector{Int}
    served::BitVector
end
get_vkey(p::YesLoadHomChargeElementaryPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.load,
    p.served,
)

mutable struct YesLoadHomChargeNoServicePPathLabel <: YesLoadHomChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNoServiceSubpathLabel}
    charging_actions::Vector{Int}
end
get_vkey(p::YesLoadHomChargeNoServicePPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.load,
)

abstract type YesLoadHomChargeNgRoutePPathLabel <: YesLoadHomChargePPathLabel end

mutable struct YesLoadHomChargeNgRouteNoCutsPPathLabel <: YesLoadHomChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNgRouteNoCutsSubpathLabel}
    charging_actions::Vector{Int}
    ng_fset::BitVector
end
get_vkey(p::YesLoadHomChargeNgRouteNoCutsPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.load,
    p.ng_fset,
)

mutable struct YesLoadHomChargeNgRouteSR3CutsPPathLabel <: YesLoadHomChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNgRouteSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::YesLoadHomChargeNgRouteSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    p.load,
    p.ng_fset,
)

mutable struct YesLoadHomChargeNgRouteLmSR3CutsPPathLabel <: YesLoadHomChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNgRouteLmSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::YesLoadHomChargeNgRouteLmSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    p.load,
    p.ng_fset,
)

abstract type YesLoadHetChargePPathLabel <: YesLoadPPathLabel end

mutable struct YesLoadHetChargeElementaryPPathLabel <: YesLoadHetChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadElementarySubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    served::BitVector
end
get_vkey(p::YesLoadHetChargeElementaryPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.load,
    - p.charge_rebalance_slacks,
    p.served,
)

mutable struct YesLoadHetChargeNoServicePPathLabel <: YesLoadHetChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNoServiceSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
end
get_vkey(p::YesLoadHetChargeNoServicePPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.load,
    - p.charge_rebalance_slacks,
)

abstract type YesLoadHetChargeNgRoutePPathLabel <: YesLoadHetChargePPathLabel end
mutable struct YesLoadHetChargeNgRouteNoCutsPPathLabel <: YesLoadHetChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNgRouteNoCutsSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    ng_fset::BitVector
end
get_vkey(p::YesLoadHetChargeNgRouteNoCutsPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.load,
    - p.charge_rebalance_slacks,
    p.ng_fset,
)

mutable struct YesLoadHetChargeNgRouteSR3CutsPPathLabel <: YesLoadHetChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNgRouteSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::YesLoadHetChargeNgRouteSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    p.load,
    - p.charge_rebalance_slacks,
    p.ng_fset,
)

mutable struct YesLoadHetChargeNgRouteLmSR3CutsPPathLabel <: YesLoadHetChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    load::Int
    subpath_labels::Vector{YesLoadNgRouteLmSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::YesLoadHetChargeNgRouteLmSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    p.load,
    - p.charge_rebalance_slacks,
    p.ng_fset,
)

abstract type NoLoadPPathLabel <: PPathLabel end

abstract type NoLoadHomChargePPathLabel <: NoLoadPPathLabel end

mutable struct NoLoadHomChargeElementaryPPathLabel <: NoLoadHomChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadElementarySubpathLabel} # TODO: replace with dict keys to nondominated subpaths
    charging_actions::Vector{Int}
    served::BitVector
end
get_vkey(p::NoLoadHomChargeElementaryPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    BitVector(p.served),
)

mutable struct NoLoadHomChargeNoServicePPathLabel <: NoLoadHomChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNoServiceSubpathLabel}
    charging_actions::Vector{Int}
end
get_vkey(p::NoLoadHomChargeNoServicePPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
)

abstract type NoLoadHomChargeNgRoutePPathLabel <: NoLoadHomChargePPathLabel end

mutable struct NoLoadHomChargeNgRouteNoCutsPPathLabel <: NoLoadHomChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNgRouteNoCutsSubpathLabel}
    charging_actions::Vector{Int}
    ng_fset::BitVector
end
get_vkey(p::NoLoadHomChargeNgRouteNoCutsPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    p.ng_fset,
)

mutable struct NoLoadHomChargeNgRouteSR3CutsPPathLabel <: NoLoadHomChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNgRouteSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::NoLoadHomChargeNgRouteSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    p.ng_fset,
)

mutable struct NoLoadHomChargeNgRouteLmSR3CutsPPathLabel <: NoLoadHomChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNgRouteLmSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::NoLoadHomChargeNgRouteLmSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    p.ng_fset,
)

abstract type NoLoadHetChargePPathLabel <: NoLoadPPathLabel end

mutable struct NoLoadHetChargeElementaryPPathLabel <: NoLoadHetChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadElementarySubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    served::BitVector
end
get_vkey(p::NoLoadHetChargeElementaryPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    - p.charge_rebalance_slacks,
    BitVector(p.served),
)

mutable struct NoLoadHetChargeNoServicePPathLabel <: NoLoadHetChargePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNoServiceSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
end
get_vkey(p::NoLoadHetChargeNoServicePPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    - p.charge_rebalance_slacks,
)

abstract type NoLoadHetChargeNgRoutePPathLabel <: NoLoadHetChargePPathLabel end
mutable struct NoLoadHetChargeNgRouteNoCutsPPathLabel <: NoLoadHetChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNgRouteNoCutsSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    ng_fset::BitVector
end
get_vkey(p::NoLoadHetChargeNgRouteNoCutsPPathLabel) = (
    p.cost,
    p.time,
    - p.charge,
    - p.charge_rebalance_slacks,
    p.ng_fset,
)

mutable struct NoLoadHetChargeNgRouteSR3CutsPPathLabel <: NoLoadHetChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNgRouteSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::NoLoadHetChargeNgRouteSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    - p.charge_rebalance_slacks,
    p.ng_fset,
)

mutable struct NoLoadHetChargeNgRouteLmSR3CutsPPathLabel <: NoLoadHetChargeNgRoutePPathLabel
    cost::Float64
    time::Int
    charge::Int
    subpath_labels::Vector{NoLoadNgRouteLmSR3CutsSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    charge_rebalance_slacks::Vector{Int}
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::NoLoadHetChargeNgRouteLmSR3CutsPPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time,
    - p.charge,
    - p.charge_rebalance_slacks,
    p.ng_fset,
)


function get_path_label_type(
    use_load::Load,
    use_time_windows::TimeWindows,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
)
    if use_load isa YesLoad
        if charge_costs isa HomCharge
            if customer_service isa Elementary
                return YesLoadHomChargeElementaryPPathLabel
            elseif customer_service isa NoService
                return YesLoadHomChargeNoServicePPathLabel
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return YesLoadHomChargeNgRouteNoCutsPPathLabel
                elseif cuts isa SR3Cuts
                    return YesLoadHomChargeNgRouteSR3CutsPPathLabel
                elseif cuts isa LmSR3Cuts
                    return YesLoadHomChargeNgRouteLmSR3CutsPPathLabel
                end
            end
        elseif charge_costs isa HetCharge
            if customer_service isa Elementary
                return YesLoadHetChargeElementaryPPathLabel
            elseif customer_service isa NoService
                return YesLoadHetChargeNoServicePPathLabel
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return YesLoadHetChargeNgRouteNoCutsPPathLabel
                elseif cuts isa SR3Cuts
                    return YesLoadHetChargeNgRouteSR3CutsPPathLabel
                elseif cuts isa LmSR3Cuts
                    return YesLoadHetChargeNgRouteLmSR3CutsPPathLabel
                end
            end
        end
    elseif use_load isa NoLoad
        if charge_costs isa HomCharge
            if customer_service isa Elementary
                return NoLoadHomChargeElementaryPPathLabel
            elseif customer_service isa NoService
                return NoLoadHomChargeNoServicePPathLabel
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return NoLoadHomChargeNgRouteNoCutsPPathLabel
                elseif cuts isa SR3Cuts
                    return NoLoadHomChargeNgRouteSR3CutsPPathLabel
                elseif cuts isa LmSR3Cuts
                    return NoLoadHomChargeNgRouteLmSR3CutsPPathLabel
                end
            end
        elseif charge_costs isa HetCharge
            if customer_service isa Elementary
                return NoLoadHetChargeElementaryPPathLabel
            elseif customer_service isa NoService
                return NoLoadHetChargeNoServicePPathLabel
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return NoLoadHetChargeNgRouteNoCutsPPathLabel
                elseif cuts isa SR3Cuts
                    return NoLoadHetChargeNgRouteSR3CutsPPathLabel
                elseif cuts isa LmSR3Cuts
                    return NoLoadHetChargeNgRouteLmSR3CutsPPathLabel
                end
            end
        end
    end
end


function create_new_subpath_label(
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
    starting_node::Int,
    data::EVRPData,
    ;
    n_cuts::Int = 0,
    neighborhoods::BitMatrix = falses(data.n_nodes, data.n_nodes),
    λmemory::BitMatrix = falses(n_cuts, data.n_nodes),
)

    if customer_service isa NgRoute
        # Π: Forward NG-set
        ng_fset = falses(data.n_nodes)
        ng_fset[starting_node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        ng_residue = copy(neighborhoods[:, starting_node])
        # Φ: Backward NG-set
        ng_bset = falses(data.n_nodes)
        ng_bset[starting_node] = true
    end
    if use_load isa YesLoad
        if customer_service isa Elementary
            return YesLoadElementarySubpathLabel(
                0.0, 0, 0, 0,
                [starting_node,], 
                falses(data.n_customers),
            )
        elseif customer_service isa NoService
            return YesLoadNoServiceSubpathLabel(
                0.0, 0, 0, 0,
                [starting_node,],
            )
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return YesLoadNgRouteNoCutsSubpathLabel(
                    0.0, 0, 0, 0,
                    [starting_node,], 
                    ng_fset, ng_residue, ng_bset,
                )
            elseif cuts isa SR3Cuts
                return YesLoadNgRouteSR3CutsSubpathLabel(
                    0.0, 0, 0, 0,
                    [starting_node,], 
                    ng_fset, ng_residue, ng_bset,
                    falses(n_cuts),
                )
            elseif cuts isa LmSR3Cuts
                return YesLoadNgRouteLmSR3CutsSubpathLabel(
                    0.0, 0, 0, 0,
                    [starting_node,], 
                    ng_fset, ng_residue, ng_bset,
                    falses(n_cuts),
                    falses(n_cuts),
                    copy(λmemory[:, starting_node]),
                )
            end
        end
    else
        if customer_service isa Elementary
            return NoLoadElementarySubpathLabel(
                0.0, 0, 0,
                [starting_node,],
                falses(data.n_customers),
            )
        elseif customer_service isa NoService
            return NoLoadNoServiceSubpathLabel(
                0.0, 0, 0,
                [starting_node,],
            )
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return NoLoadNgRouteNoCutsSubpathLabel(
                    0.0, 0, 0,
                    [starting_node,], 
                    ng_fset, ng_residue, ng_bset,
                )
            elseif cuts isa SR3Cuts
                return NoLoadNgRouteSR3CutsSubpathLabel(
                    0.0, 0, 0,
                    [starting_node,], 
                    ng_fset, ng_residue, ng_bset,
                    falses(n_cuts),
                )
            elseif cuts isa LmSR3Cuts
                return NoLoadNgRouteLmSR3CutsSubpathLabel(
                    0.0, 0, 0,
                    [starting_node,], 
                    ng_fset, ng_residue, ng_bset,
                    falses(n_cuts),
                    falses(n_cuts),
                    copy(λmemory[:, starting_node]),
                )
            end
        end
    end
end

function create_empty_path_label(
    use_load::Load,
    use_time_windows::TimeWindows,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
    depot::Int,
    data::EVRPData,
    ;
    n_cuts::Int = 0,
    neighborhoods::BitMatrix = falses(data.n_nodes, data.n_nodes),
    λmemory::BitMatrix = falses(n_cuts, data.n_nodes),
)
    s = create_new_subpath_label(
        use_load,
        use_time_windows,
        customer_service,
        cuts,
        depot,
        data,
        ;
        n_cuts = n_cuts,
        neighborhoods = neighborhoods,
        λmemory = λmemory,
    )
    push!(s.nodes, depot)
    p = create_new_path_label(
        use_load,
        charge_costs,
        customer_service,
        cuts,
        depot,
        data,
        ;
        n_cuts = n_cuts,
    )
    push!(p.subpath_labels, s)

    return p
end

function update_subpath_ngroute!(
    new_subpath::SubpathLabel,
    next_node::Int,
    neighborhoods::BitMatrix
)
    # update forward ng-set
    ngroute_update_fset!(new_subpath, next_node, neighborhoods)
    ## IMPORTANT: update backward ng-set before residue, 
    ## since residue affects backward ng-set
    # update backward ng-set
    ngroute_update_bset!(new_subpath, next_node)
    # update ng-residue
    ngroute_update_residue!(new_subpath, next_node, neighborhoods)
    return 
end

function update_subpath_cut_labels_SR3!(
    new_subpath::SubpathLabel,
    next_node::Int,
    λvals::Vector{Float64},
    λcust::BitMatrix,
)
    ## IMPORTANT: update cost before flabels
    ## since flabels affect cost
    SR3_update_cost!(new_subpath, next_node, λvals, λcust)
    SR3_update_cut_flabels!(new_subpath, next_node, λcust)
    return
end

function update_subpath_cut_labels_LmSR3!(
    new_subpath::SubpathLabel,
    next_node::Int,
    λvals::Vector{Float64},
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    ## IMPORTANT: update cost before flabels
    ## since flabels affect cost
    lmSR3_update_cost!(new_subpath, next_node, λvals, λcust, λmemory)
    # Update cut flabels
    lmSR3_update_cut_flabels!(new_subpath, next_node, λcust, λmemory)
    # Update cut blabels
    lmSR3_update_cut_blabels!(new_subpath, next_node, λcust)
    # Update cut mlabels
    lmSR3_update_cut_mlabels!(new_subpath, next_node, λmemory)
    return
end

function compute_new_subpath(
    current_subpath::SubpathLabel,
    graph::EVRPGraph,
    modified_costs::Matrix{Float64},
    current_node::Int,
    next_node::Int,
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
    ;
    neighborhoods::BitMatrix = falses(graph.n_nodes, graph.n_nodes),
    n_cuts::Int = 0,
    λvals::Vector{Float64} = Float64[],
    λcust::BitMatrix = falses(n_cuts, graph.n_nodes),
    λmemory::BitMatrix = falses(n_cuts, graph.n_nodes),
)
    # customer service
    if customer_service isa Elementary
        if next_node in graph.N_customers && current_subpath.served[next_node]
            return (false, current_subpath)
        end
    elseif customer_service isa NgRoute
        if current_subpath.ng_fset[next_node]
            return (false, current_subpath)
        end
    end

    # time and charge feasibility
    if current_subpath.time_taken + graph.t[current_node, next_node] + graph.min_t[next_node] > graph.T
        return (false, current_subpath)
    end

    if current_subpath.charge_taken + graph.q[current_node, next_node] + graph.min_q[next_node] > graph.B
        return (false, current_subpath)
    end

    # load feasibility
    if use_load isa YesLoad
        if current_subpath.load + graph.d[next_node] > graph.C
            return (false, current_subpath)
        end
    end

    # Create new subpath
    new_subpath = copy(current_subpath)

    new_subpath.cost += modified_costs[current_node, next_node]
    new_subpath.time_taken += graph.t[current_node, next_node]
    new_subpath.charge_taken += graph.q[current_node, next_node]
    push!(new_subpath.nodes, next_node)

    # Load
    if use_load isa YesLoad
        new_subpath.load += graph.d[next_node]
    end

    # Customer service
    if customer_service isa Elementary
        if next_node in graph.N_customers
            new_subpath.served[next_node] += 1
        end
    elseif customer_service isa NgRoute
        update_subpath_ngroute!(new_subpath, next_node, neighborhoods)
    end

    # Cuts
    if cuts isa SR3Cuts
        update_subpath_cut_labels_SR3!(new_subpath, next_node, λvals, λcust)
    elseif cuts isa LmSR3Cuts 
        update_subpath_cut_labels_LmSR3!(new_subpath, next_node, λvals, λcust, λmemory)
    end

    return (true, new_subpath)
end


function generate_subpath_labels(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
    ;
    neighborhoods::BitMatrix = falses(graph.n_nodes, graph.n_nodes),
    λ::Dict{<:Union{
        NTuple{3, Int},
        Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
    }, Float64} = Dict{NTuple{3, Int}, Float64}(),
    time_limit::Float64 = Inf,
)
    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    if cuts isa NoCuts
        λvals = Float64[]
        λcust = falses(length(λ), graph.n_nodes)
        λmemory = falses(length(λ), graph.n_nodes)
    elseif cuts isa SR3Cuts
        λvals, λcust = prepare_lambda(λ, graph.n_nodes)
        λmemory = falses(length(λ), graph.n_nodes)
    elseif cuts isa LmSR3Cuts
        λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)
    end

    subpath_labels = Dict(
        (starting_node, current_node) => SortedDict{
            get_subpath_vkey_type(use_load, customer_service, cuts),
            get_subpath_label_type(use_load, customer_service, cuts),
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )
    unexplored_states = SortedSet{get_subpath_vkey_fkey_type(use_load, customer_service, cuts)}()

    for starting_node in graph.N_depots_charging
        s = create_new_subpath_label(
            use_load,
            use_time_windows,
            customer_service,
            cuts,
            starting_node,
            data,
            ;
            n_cuts = length(λ), 
            neighborhoods = neighborhoods,
            λmemory = λmemory, 
        )
        vkey = get_vkey(s)
        subpath_labels[(starting_node, starting_node)][vkey] = s
        push!(unexplored_states, (vkey..., starting_node, starting_node,))
    end
    while length(unexplored_states) > 0
        # Time limit check
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        # Retrieve most promising unexplored state
        state = pop!(unexplored_states)
        (current_vkey..., starting_node, current_node) = state
        # If current vkey is already stale, skip
        if !(current_vkey in keys(subpath_labels[(starting_node, current_node)]))
            continue
        end

        current_subpath = subpath_labels[(starting_node, current_node)][current_vkey]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            # Check feasibility and create new subpath
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath,
                graph,
                modified_costs,
                current_node,
                next_node, 
                use_load,
                use_time_windows,
                customer_service,
                cuts,
                ;
                neighborhoods = neighborhoods,
                n_cuts = length(λ),
                λvals = λvals,
                λcust = λcust,
                λmemory = λmemory,
            )
            !feasible && continue
            
            # Create keys, add to collection, and add to unexplored states
            new_vkey = get_vkey(new_subpath)
            if cuts isa NoCuts
                added = add_label_to_collection!(
                    subpath_labels[(starting_node, next_node)],
                    new_vkey, new_subpath,
                )
            elseif cuts isa SR3Cuts
                added = add_label_to_collection_cuts!(
                    subpath_labels[(starting_node, next_node)],
                    new_vkey, new_subpath, λvals,
                )
            elseif cuts isa LmSR3Cuts
                added = add_label_to_collection_lmSR3_subpath!(
                    subpath_labels[(starting_node, next_node)],
                    new_vkey, new_subpath, λvals,
                )
            end
            if added && next_node in graph.N_customers
                push!(unexplored_states, (new_vkey..., starting_node, next_node,))
            end
        end
    end

    # Cleanup
    for starting_node in graph.N_depots_charging, end_node in graph.N_customers
        pop!(subpath_labels, (starting_node, end_node))
    end
    for node in graph.N_depots_charging
        s = create_new_subpath_label(
            use_load,
            use_time_windows,
            customer_service,
            cuts,
            node,
            data,
            ;
            n_cuts = length(λ),
            neighborhoods = neighborhoods,
            λmemory = λmemory,
        )
        pop!(subpath_labels[(node, node)], get_vkey(s))
    end
    for starting_node in graph.N_depots, end_node in graph.N_depots_charging
        for s in values(subpath_labels[(starting_node, end_node)])
            s.cost -= κ[starting_node]
        end
    end
    for starting_node in graph.N_depots_charging, end_node in graph.N_depots
        for s in values(subpath_labels[(starting_node, end_node)])
            s.cost -= μ[end_node]
        end
    end

    return subpath_labels

end


function create_new_path_label(
    use_load::Load,
    use_time_windows::TimeWindows,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
    starting_node::Int,
    data::EVRPData,
    ;
    n_cuts::Int = 0,
)
    if customer_service isa NgRoute
        ng_fset = falses(data.n_nodes)
        ng_fset[starting_node] = true
    end

    if use_load isa YesLoad
        if charge_costs isa HomCharge
            if customer_service isa Elementary
                return YesLoadHomChargeElementaryPPathLabel(
                    0.0, 0, data.B, 0,
                    YesLoadElementarySubpathLabel[], 
                    Int[],
                    falses(data.n_customers),
                )
            elseif customer_service isa NoService
                return YesLoadHomChargeNoServicePPathLabel(
                    0.0, 0, data.B, 0,
                    YesLoadNoServiceSubpathLabel[], 
                    Int[],
                )
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return YesLoadHomChargeNgRouteNoCutsPPathLabel(
                        0.0, 0, data.B, 0,
                        YesLoadNgRouteNoCutsSubpathLabel[], 
                        Int[],
                        ng_fset,
                    )
                elseif cuts isa SR3Cuts
                    return YesLoadHomChargeNgRouteSR3CutsPPathLabel(
                        0.0, 0, data.B, 0,
                        YesLoadNgRouteSR3CutsSubpathLabel[], 
                        Int[],
                        ng_fset, falses(n_cuts),
                    )
                elseif cuts isa LmSR3Cuts
                    return YesLoadHomChargeNgRouteLmSR3CutsPPathLabel(
                        0.0, 0, data.B, 0,
                        YesLoadNgRouteLmSR3CutsSubpathLabel[], 
                        Int[],
                        ng_fset, falses(n_cuts),
                    )
                end
            end
        elseif charge_costs isa HetCharge
            if customer_service isa Elementary
                return YesLoadHetChargeElementaryPPathLabel(
                    0.0, 0, data.B, 0,
                    YesLoadElementarySubpathLabel[], 
                    Int[], 
                    zeros(Int, data.charge_cost_nlevels), 
                    zeros(Int, data.charge_cost_nlevels - 1),
                    falses(data.n_customers),
                )
            elseif customer_service isa NoService
                return YesLoadHetChargeNoServicePPathLabel(
                    0.0, 0, data.B, 0,
                    YesLoadNoServiceSubpathLabel[], 
                    Int[], 
                    zeros(Int, data.charge_cost_nlevels), 
                    zeros(Int, data.charge_cost_nlevels - 1),
                )
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return YesLoadHetChargeNgRouteNoCutsPPathLabel(
                        0.0, 0, data.B, 0,
                        YesLoadNgRouteNoCutsSubpathLabel[], 
                        Int[], 
                        zeros(Int, data.charge_cost_nlevels), 
                        zeros(Int, data.charge_cost_nlevels - 1),
                        ng_fset, 
                    )
                elseif cuts isa SR3Cuts
                    return YesLoadHetChargeNgRouteSR3CutsPPathLabel(
                        0.0, 0, data.B, 0,
                        YesLoadNgRouteSR3CutsSubpathLabel[], 
                        Int[], 
                        zeros(Int, data.charge_cost_nlevels), 
                        zeros(Int, data.charge_cost_nlevels - 1),
                        ng_fset, falses(n_cuts),
                    )
                elseif cuts isa LmSR3Cuts
                    return YesLoadHetChargeNgRouteLmSR3CutsPPathLabel(
                        0.0, 0, data.B, 0,
                        YesLoadNgRouteLmSR3CutsSubpathLabel[], 
                        Int[], 
                        zeros(Int, data.charge_cost_nlevels), 
                        zeros(Int, data.charge_cost_nlevels - 1),
                        ng_fset, falses(n_cuts),
                    )
                end
            end
        end
    elseif use_load isa NoLoad
        if charge_costs isa HomCharge
            if customer_service isa Elementary
                return NoLoadHomChargeElementaryPPathLabel(
                    0.0, 0, data.B,
                    NoLoadElementarySubpathLabel[], 
                    Int[],
                    falses(data.n_customers),
                )
            elseif customer_service isa NoService
                return NoLoadHomChargeNoServicePPathLabel(
                    0.0, 0, data.B,
                    NoLoadNoServiceSubpathLabel[], 
                    Int[],
                )
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return NoLoadHomChargeNgRouteNoCutsPPathLabel(
                        0.0, 0, data.B,
                        NoLoadNgRouteNoCutsSubpathLabel[], 
                        Int[],
                        ng_fset,
                    )
                elseif cuts isa SR3Cuts
                    return NoLoadHomChargeNgRouteSR3CutsPPathLabel(
                        0.0, 0, data.B,
                        NoLoadNgRouteSR3CutsSubpathLabel[], 
                        Int[],
                        ng_fset, falses(n_cuts),
                    )
                elseif cuts isa LmSR3Cuts
                    return NoLoadHomChargeNgRouteLmSR3CutsPPathLabel(
                        0.0, 0, data.B,
                        NoLoadNgRouteLmSR3CutsSubpathLabel[], 
                        Int[],
                        ng_fset, falses(n_cuts),
                    )
                end
            end
        elseif charge_costs isa HetCharge
            if customer_service isa Elementary
                return NoLoadHetChargeElementaryPPathLabel(
                    0.0, 0, data.B,
                    NoLoadElementarySubpathLabel[], 
                    Int[], 
                    zeros(Int, data.charge_cost_nlevels), 
                    zeros(Int, data.charge_cost_nlevels - 1),
                    falses(data.n_customers),
                )
            elseif customer_service isa NoService
                return NoLoadHetChargeNoServicePPathLabel(
                    0.0, 0, data.B,
                    NoLoadNoServiceSubpathLabel[], 
                    Int[], 
                    zeros(Int, data.charge_cost_nlevels), 
                    zeros(Int, data.charge_cost_nlevels - 1),
                )
            elseif customer_service isa NgRoute
                if cuts isa NoCuts
                    return NoLoadHetChargeNgRouteNoCutsPPathLabel(
                        0.0, 0, data.B,
                        NoLoadNgRouteNoCutsSubpathLabel[], 
                        Int[], 
                        zeros(Int, data.charge_cost_nlevels), 
                        zeros(Int, data.charge_cost_nlevels - 1),
                        ng_fset, 
                    )
                elseif cuts isa SR3Cuts
                    return NoLoadHetChargeNgRouteSR3CutsPPathLabel(
                        0.0, 0, data.B,
                        NoLoadNgRouteSR3CutsSubpathLabel[], 
                        Int[], 
                        zeros(Int, data.charge_cost_nlevels), 
                        zeros(Int, data.charge_cost_nlevels - 1),
                        ng_fset, falses(n_cuts),
                    )
                elseif cuts isa LmSR3Cuts
                    return NoLoadHetChargeNgRouteLmSR3CutsPPathLabel(
                        0.0, 0, 0,
                        NoLoadNgRouteLmSR3CutsSubpathLabel[], 
                        Int[], 
                        zeros(Int, data.charge_cost_nlevels), 
                        zeros(Int, data.charge_cost_nlevels - 1),
                        ng_fset, falses(n_cuts),
                    )
                end
            end
        end
    end
end

function check_path_ngroute(
    path_fset::BitVector,
    node::Int,
    subpath_fset::BitVector,
    subpath_residue::BitVector,
    subpath_bset::BitVector,
)
    new_path_inf = path_fset .& subpath_bset
    new_path_inf[node] = false
    if any(new_path_inf)
        return (false, path_fset)
    else
        return (true, subpath_fset .| (path_fset .& subpath_residue))
    end
end

function update_path_cut_labels_SR3!(
    new_path::PPathLabel,
    s::SubpathLabel,
    λvals::Vector{Float64},
)
    # Update reduced cost
    new_path.cost -= sum(λvals[new_path.cut_flabels .& s.cut_flabels])
    # Update cut labels
    new_path.cut_flabels .⊻= s.cut_flabels
    return
end

function update_path_cut_labels_lmSR3!(
    new_path::PPathLabel,
    s::SubpathLabel,
    λvals::Vector{Float64},
)
    # Update reduced cost
    new_path.cost -= sum(λvals[new_path.cut_flabels .& s.cut_blabels])
    # Update cut labels
    new_path.cut_flabels .&= s.cut_mlabels
    new_path.cut_flabels .⊻= s.cut_flabels
    return
end

function compute_new_path(
    current_path::PPathLabel,
    data::EVRPData,
    graph::EVRPGraph,
    current_node::Int,
    next_node::Int,
    s::SubpathLabel,
    use_load::Load,
    use_time_windows::TimeWindows,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
    ;
    λvals::Vector{Float64} = Float64[],
)

    # don't count initial subpath again
    if (
        next_node in graph.N_depots
        && s.time_taken == 0
    )
        return (false, current_path)
    end

    # elementarity
    if customer_service isa Elementary
        if any(s.served + current_path.served .> 1)
            return (false, current_path)
        end
    elseif customer_service isa NgRoute
        (feasible, new_path_ng_fset) = check_path_ngroute(
            current_path.ng_fset, 
            current_node,
            s.ng_fset,
            s.ng_residue,
            s.ng_bset,
        )
        if !feasible
            return (false, current_path)
        end
    end

    # time and charge feasibility
    (original_charge_amount, end_time, end_charge) = charge_to_specified_level(
        current_path.charge,
        s.charge_taken, # desired charge  
        current_path.time,
    )
    end_time += s.time_taken
    end_charge -= s.charge_taken
    if end_time + graph.min_t[next_node] > graph.T
        return (false, current_path)
    end

    # load feasibility
    if use_load isa YesLoad
        end_load = current_path.load + s.load
        if end_load > graph.C
            return (false, current_path)
        end
    end

    # Create new path
    new_path = copy(current_path)
    new_path.cost += s.cost
    new_path.time = end_time
    new_path.charge = end_charge
    if length(current_path.subpath_labels) > 0
        push!(new_path.charging_actions, original_charge_amount)
        new_path.cost += data.charge_cost_coeff * original_charge_amount
    end

    if use_load isa YesLoad
        new_path.load = end_load
    end

    push!(new_path.subpath_labels, s)

    # Customer service
    if customer_service isa Elementary
        new_path.served .+= s.served
    elseif customer_service isa NgRoute
        new_path.ng_fset = new_path_ng_fset
    end

    # Cuts
    if cuts isa SR3Cuts
        update_path_cut_labels_SR3!(new_path, s, λvals)
    elseif cuts isa LmSR3Cuts
        update_path_cut_labels_lmSR3!(new_path, s, λvals)
    end

    return (true, new_path)    
end


function generate_path_labels(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    use_load::Load,
    use_time_windows::TimeWindows,
    charge_costs::ChargeCosts,
    customer_service::CustomerService,
    cuts::Cuts,
    subpath_labels::Dict{
        Tuple{Int, Int},
        <:AbstractDict
    },
    ;
    λ::Dict{<:Tuple, Float64} = Dict(),
    time_limit::Float64 = Inf,
)
    start_time = time()
    if cuts isa NoCuts
        λvals = Float64[]
        λmemory = falses(length(λ), graph.n_nodes)
    end
    if cuts isa SR3Cuts
        λvals, _ = prepare_lambda(λ, graph.n_nodes)
        λmemory = falses(length(λ), graph.n_nodes)
    elseif cuts isa LmSR3Cuts
        λvals, _, λmemory = prepare_lambda(λ, graph.n_nodes)
    end

    path_labels = Dict(
        (starting_node, current_node) => SortedDict{
            get_path_vkey_type(use_load, use_time_windows, charge_costs, customer_service, cuts),
            get_path_label_type(use_load, use_time_windows, charge_costs, customer_service, cuts),
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )
    unexplored_states = SortedSet{get_path_vkey_fkey_type(use_load, use_time_windows, charge_costs, customer_service, cuts)}()

    for depot in graph.N_depots
        p = create_new_path_label(
            use_load,
            use_time_windows,
            charge_costs,
            customer_service,
            cuts,
            depot,
            data,
            ;
            n_cuts = length(λ),
        )
        vkey = get_vkey(p)
        path_labels[(depot, depot)][vkey] = p
        push!(unexplored_states, (vkey..., depot, depot,))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        # Retrieve most promising unexplored state
        state = pop!(unexplored_states)
        (current_vkey..., starting_node, current_node) = state
        if !(current_vkey in keys(path_labels[(starting_node, current_node)]))
            continue
        end
        current_path = path_labels[(starting_node, current_node)][current_vkey]
        for next_node in graph.N_depots_charging
            for (_, s) in pairs(subpath_labels[current_node, next_node])
                (feasible, new_path) = compute_new_path(
                    current_path,
                    data,
                    graph,
                    current_node,
                    next_node,
                    s,
                    use_load,
                    use_time_windows,
                    charge_costs,
                    customer_service,
                    cuts,
                    ;
                    λvals = λvals,
                )
                !feasible && continue
                
                new_vkey = get_vkey(new_path)
                if cuts isa NoCuts
                    added = add_label_to_collection!(
                        path_labels[(starting_node, next_node)],
                        new_vkey, new_path,
                    )
                elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                    added = add_label_to_collection_cuts!(
                        path_labels[(starting_node, next_node)],
                        new_vkey, new_path, λvals,
                    )
                end
                if added && next_node in graph.N_charging
                    push!(unexplored_states, (new_vkey..., starting_node, next_node,))
                end
            end
        end
    end

    for depot in graph.N_depots
        p = create_empty_path_label(
            use_load,
            use_time_windows,
            charge_costs,
            customer_service,
            cuts,
            depot,
            data,
            ;
            n_cuts = length(λ),
            λmemory = λmemory, 
        )
        vkey = get_vkey(p)
        p.subpath_labels[1].cost -= (κ[depot] + μ[depot])
        p.cost -= (κ[depot] + μ[depot])
        path_labels[(depot, depot)][vkey] = p
    end

    for starting_node in graph.N_depots, end_node in graph.N_charging
        pop!(path_labels, (starting_node, end_node))
    end

    return path_labels

end


unwrap_path_labels(p::Label) = Label[p]

function unwrap_path_labels(d::AbstractDict)
    u = Label[]
    for v in values(d)
        append!(u, unwrap_path_labels(v))
    end
    return u
end

function get_negative_path_labels_from_path_labels(
    path_labels::Dict{
        Tuple{Int, Int}, 
        T,
    },
) where {T <: AbstractDict}
    return Label[
        path_label
        for path_label in unwrap_path_labels(path_labels)
            if path_label.cost < -1e-4
    ]
end

# function subtract_duals!(
#     path_labels::Dict{
#         Tuple{Int, Int}, 
#         <:AbstractDict,
#     },
#     κ::Dict{Int, Float64},
#     μ::Dict{Int, Float64},
# )
#     for (starting_node, end_node) in keys(path_labels)
#         # println(starting_node, " -> ", end_node)
#         # println([p.cost for p in values(path_labels[(starting_node, end_node)])])
#         for (vkey, p) in pairs(path_labels[(starting_node, end_node)])
#             p.cost -= (κ[starting_node] + μ[end_node])
#         end
#         # println([p.cost for p in values(path_labels[(starting_node, end_node)])])
#     end
#     return path_labels
#     return
# end



function subproblem_iteration_ours(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{<:Any, Float64},
    ;
    load::Bool = false,
    time_windows::Bool = false,
    charge_cost_heterogenous::Bool = false,
    neighborhoods::BitMatrix = falses(graph.n_nodes, graph.n_nodes),
    ngroute::Bool = false,
    elementary::Bool = true,
    time_limit::Float64 = Inf,
)
    start_time = time()

    use_load = load ? YesLoad() : NoLoad()
    use_time_windows = time_windows ? YesTimeWindows() : NoTimeWindows()
    charge_costs = charge_cost_heterogenous ? HetCharge() : HomCharge()
    if ngroute
        customer_service = NgRoute()
        if length(λ) == 0
            cuts = NoCuts()
        elseif keytype(λ) == NTuple{3, Int}
            cuts = SR3Cuts()
        elseif keytype(λ) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
            cuts = LmSR3Cuts()
        else
            error("Unrecognized key type for λ: $(keytype(λ))")
        end
    else
        neighborhoods = falses(graph.n_nodes, graph.n_nodes)
        customer_service = elementary ? Elementary() : NoService()
        cuts = NoCuts()
    end

    subpath_labels_result = @timed generate_subpath_labels(
        data,
        graph,
        κ,
        μ,
        ν,
        use_load,
        use_time_windows,
        customer_service,
        cuts,
        ;
        neighborhoods = neighborhoods,
        λ = λ,
        time_limit = time_limit - (time() - start_time),
    )
    subpath_labels_time = subpath_labels_result.time

    path_labels_result = @timed generate_path_labels(
        data,
        graph,
        κ,
        μ,
        use_load,
        use_time_windows,
        charge_costs,
        customer_service,
        cuts,
        subpath_labels_result.value,
        ;
        λ = λ,
        time_limit = time_limit - (time() - start_time),
    )
    path_labels_time = path_labels_result.time
    path_labels = path_labels_result.value
    # subtract_duals!(path_labels, κ, μ)

    negative_path_labels = get_negative_path_labels_from_path_labels(path_labels)
    negative_path_labels_count = length(negative_path_labels)

    return (
        negative_path_labels,
        negative_path_labels_count,
        subpath_labels_time,
        path_labels_time,
    )
end

function convert_path_label_to_path(
    path_label::PPathLabel,
    data::EVRPData,
    graph::EVRPGraph,
    ;
    use_load::Bool = false,
)
    current_time, current_charge = (0, graph.B)
    prev_time, prev_charge = current_time, current_charge
    s_labels = copy(path_label.subpath_labels)
    deltas = copy(path_label.charging_actions)
    p = Path(
        subpaths = Subpath[],
        charging_arcs = ChargingArc[],
        served = zeros(Int, graph.n_customers),
        load = 0,
        arcs = NTuple{2, Int}[],
        customer_arcs = NTuple{2, Int}[],
    )
    while true
        s_label = popfirst!(s_labels)
        prev_time = current_time
        prev_charge = current_charge
        current_node = s_label.nodes[end]
        current_time = current_time + s_label.time_taken
        current_charge = current_charge - s_label.charge_taken
        served = [count(x -> x == i, s_label.nodes[2:end-1]) for i in 1:graph.n_customers]
        s = Subpath(
            n_customers = graph.n_customers,
            starting_node = s_label.nodes[1],
            starting_time = prev_time,
            starting_charge = prev_charge,
            current_node = current_node,
            arcs = collect(zip(s_label.nodes[1:end-1], s_label.nodes[2:end])),
            current_time = current_time,
            current_charge = current_charge,
            load = (use_load ? s_label.load : 0),
            served = served,
        )
        push!(p.subpaths, s)
        if length(deltas) == 0 
            break
        end
        delta = popfirst!(deltas)
        prev_time = current_time
        prev_charge = current_charge
        current_time = current_time + delta
        current_charge = current_charge + delta
        a = ChargingArc(
            node = current_node,
            time_start = prev_time,
            time_end = current_time,
            time_diff = delta,
            charge_start = prev_charge,
            charge_end = current_charge,
            charge_diff = delta,
        )
        push!(p.charging_arcs, a)
    end
    p.served = sum(s.served for s in p.subpaths)
    p.load = sum(s.load for s in p.subpaths)
    p.arcs = vcat([s.arcs for s in p.subpaths]...)
    customers = [a[1] for a in p.arcs if a[1] in graph.N_customers]
    p.customer_arcs = collect(zip(customers[1:end-1], customers[2:end]))
    return p
end
