# the real case uses municipality instead of global setting

#=
    The  GEMS framework makes it very easy to include
    various types of transmission rates, detemining the
    likelyhood of a disease transmission dependent on the
    characteristics of the involved agents and the setting. 
=#
#using CSV
# first, load the GEMS package and required dependencies
# Note: if you define a new plot within the GEMS package,
# all these dependencies will be readily available
using GEMS, Parameters, DataFrames, TOML, Plots, FileIO, Distributions, CSV, CategoricalArrays, JLD2, Random
BASE_FOLDER = dirname(dirname(pathof(GEMS)))

##################### Distribution ###########################
using Distributions, StatsFuns

function GEMS.create_distribution(params::Vector{Float64}, type::String)
    try 
        return eval(Symbol(type))(params...)
    catch
        throw("Distribution $type cannot be created with parameters: $params")
    end
end

struct ZeroTruncatedPoisson <: DiscreteUnivariateDistribution
    λ::Float64
    function ZeroTruncatedPoisson(λ::Float64)
        λ > 0 || throw(ArgumentError("λ must be positive."))
        new(λ)
    end
end

function ZeroTruncatedPoisson(params)
    length(params) == 1 || throw(ArgumentError("Expected one parameter for ZeroTruncatedPoisson"))
    return ZeroTruncatedPoisson(float(params[1]))
end

# Support
Distributions.support(d::ZeroTruncatedPoisson) = 1:typemax(Int)

# PMF
function Distributions.pdf(d::ZeroTruncatedPoisson, k::Int)
    k < 1 && return 0.0
    λ = d.λ
    pois = Poisson(λ)
    return pdf(pois, k) / (1 - pdf(pois, 0))
end

# log PMF
function Distributions.logpdf(d::ZeroTruncatedPoisson, k::Int)
    k < 1 && return -Inf
    λ = d.λ
    pois = Poisson(λ)
    return logpdf(pois, k) - log(1 - pdf(pois, 0))
end

# Sampling
function Base.rand(d::ZeroTruncatedPoisson)
    λ = d.λ
    while true
        x = rand(Poisson(λ))
        x > 0 && return x
    end
end

# Mean
function Distributions.mean(d::ZeroTruncatedPoisson)
    λ = d.λ
    return λ / (1 - exp(-λ))
end

# Variance
function Distributions.var(d::ZeroTruncatedPoisson)
    λ = d.λ
    Z = 1 - exp(-λ)
    return λ * (1 + λ * exp(-λ) / Z) / Z - (λ / Z)^2
end

# CDF
function Distributions.cdf(d::ZeroTruncatedPoisson, k::Int)
    k < 1 && return 0.0
    λ = d.λ
    pois = Poisson(λ)
    return (cdf(pois, k) - pdf(pois, 0)) / (1 - pdf(pois, 0))
end


d = ZeroTruncatedPoisson(2.822)

mean(d)        # Expected value
#var(d)         # Variance
#cdf(d, 5)      # Cumulative probability up to 5
#rand(d, 10)    # Generate 10 random samples
###############################################################


### SETTING UP THE NEW STRUCT
# A new TransmissionFunction requires the definition of a struct
# that contains the parameters relevant for the calculation of the
# transmission rate.
# Here we define an imaginary disease that has distinct likelyhoods
# of transmission between males, between male and female and between
# females.
# Structs can be defined using the @with_kw keyword that
# automatically creates a constructor which takes the values of 
# the Sturcts parameters as kwargs. If however some additional calculations
# need to be performed in the constructor the struct should be defined
# without the @with_kw flag and a constructor must be defined by the user.

@with_kw mutable struct SettingSpecificTransmissionRate <: TransmissionFunction

    # Here we do not define a constructor but use the @with_kw flag

    # Define distribution and parameters
    distribution::String = "Normal"
    hhParameter::Vector = [0.01, 0.0001]
    wpParameter::Vector = [0.005, 0.0001]
    muParameter::Vector = [0.005, 0.0001]
    scParameter::Vector = [0.005, 0.0001]
    # Define the specific distribtions from this
    mmDistribution::Distribution = eval(Meta.parse(distribution))(hhParameter...)
    mfDistribution::Distribution = eval(Meta.parse(distribution))(wpParameter...)
    ffDistribution::Distribution = eval(Meta.parse(distribution))(muParameter...)
    scDistribution::Distribution = eval(Meta.parse(distribution))(scParameter...)

end


function immunity(tick::Int16, Ind::Individual)::Float64
    if  -1 < removed_tick(Ind) <= tick # if the agent has already recovered 100% 
        return 1.0
    else
        return 0.0
    end   
end
# Each TransmissionFunction requires the transmission_probability 
# function which is applied during the simulation to actually calculate
# the transmission rate for a specific contact. The function takes a TransmissionFunction 
# of a specific type as the first argument and thereby can be defined for different 
# Transmission functions in different ways.
# Here we determine if the contact is mm, mf or ff and return a sample from the corresponding
# distribution. 

function GEMS.transmission_probability(transFunc::SettingSpecificTransmissionRate, infecter::Individual, infected::Individual, setting::Setting, tick::Int16)::Float64
    if isa(setting, Household)
        return rand(transFunc.mmDistribution)*(1-immunity(tick, infected))
    elseif isa(setting, Office)
        return rand(transFunc.mfDistribution)*(1-immunity(tick, infected))
    elseif isa(setting, GlobalSetting)
        return rand(transFunc.ffDistribution)*(1-immunity(tick, infected))
    elseif isa(setting, SchoolClass)
        return rand(transFunc.scDistribution)*(1-immunity(tick, infected))
    elseif isa(setting, Municipality)
        return rand(transFunc.ffDistribution)*(1-immunity(tick, infected))
    else
        return 0.0
    end
end

### TESTING
# Additionaly to the creation of new transmission function appropirate tests should be defined
# in tests/infectionstest.jl. There is a specific part for testing each transmission struct and function

### RUNNING A SIMULATION AND PLOTTING RESULTS


######################################## Contacts Sampling ############################################

@with_kw struct FixedContacts <: ContactSamplingMethod
    distribution::String = "Poisson"
    mean_number_of_contacts_weekday::Float64
    mean_number_of_contacts_weekend::Float64
    #mean_number_of_contacts_friday::Float64
    wdDistribution::Distribution = eval(Meta.parse(distribution))(mean_number_of_contacts_weekday...)
    weDistribution::Distribution = eval(Meta.parse(distribution))(mean_number_of_contacts_weekend...)
    #fDistribution::Distribution = eval(Meta.parse(distribution))(mean_number_of_contacts_friday...)
    # the fact, that the contacts differ between weekdays and weekend will be implement in the `sample_contacts()` function 
end

#=
create an "override" (type in "function GEMS.sample_contacts() ... end") of the "sample_contacts()" function. This function will contain your specific sampling logic and take your created struct as an input.
it's important to note, that every custom "sample_contacts()" function needs the argumetns of type (::ContactSamplingMethod, ::Setting, ::Individual, ::Vector{Individual}, ::Int16) (the ::Int16 is the current "tick" of the simulation. When using the "contact_samples()" method, this defaults to "0") so that the simulation can automatically use the new function!

Now we define the beforementioned "sample_contacts()" function. This function incorporates the "sampling logic" for our struct "FixedContacts". Here we want to sample contacts based on the current tick and a pre-defined number of contacts. To connect the "ticks" of the Simulation to a "real world calender" we define "tick 0" as a monday.
Following this a "weekend" would occur every 6 and 7 ticks (0 = Mon, 1 = Tue, ... , 5 = Sat, 6 = Sun). We can use this information to create a function that checks, whether a tick is a "weekday" or "weekend" by taking the modulo of 7 of the tick.
=#


"""

# Parameters:
- `fixed_contacts` = your own struct
- `setting` = setting of the `ego`
- `ego` = Individual, for which the contacts should be sampled
- `present_individuals` = Individuals, currently present in `setting` 
- `tick` = current tick of the simulation object
"""
function GEMS.sample_contacts(fixed_contacts::FixedContacts, setting::Setting, ego::Individual, present_individuals::Vector{Individual},tick::Int16)
    rand(fixed_contacts.wdDistribution)
    # get the parameters stored in the "FixedContacts" struct
    # contacts_weekday::Int64 = rand(fixed_contacts.wdDistribution)
    # contacts_weekend::Int64 = rand(fixed_contacts.weDistribution)
    # contacts_friday::Int64 = rand(fixed_contacts.fDistribution)

    if isempty(present_individuals)
        throw(ArgumentError("No Individual is present in $setting. Please provide a Setting, where at least 1 Individual is present!"))
    end

    if length(present_individuals) == 1
        return Individual[]
    end

    num_of_contacts = 0

    if (tick % 7) == 5 || (tick % 7) == 6
        num_of_contacts = rand(fixed_contacts.weDistribution)
    # elseif (tick % 7) == 4 #friday
    #     num_of_contacts = rand(fixed_contacts.fDistribution)
    else
        num_of_contacts = rand(fixed_contacts.wdDistribution)
    end

    res = Vector{Individual}(undef, num_of_contacts)

    cnt = 0
    # Draw until contact list is filled, skip whenever the index individual was selected
    while cnt < num_of_contacts
        contact = rand(present_individuals)
        # if contact is NOT index individual, add them to contact list
        if Ref(contact) .!== Ref(ego)
            res[cnt + 1] = contact
            cnt += 1
        end
    end

    return res
end

#=
For our second struct we define "HouseholdSizeBasedContactSampling". This ContactSamplingMethod will sample contacts based on the size of the Household (the greater the Household size, the more contacts will be sampled). Since we can infer the Household size directly in the Simulation, we don't need to define any parameters here.
=#
struct HouseholdSizeBasedContactSampling <: ContactSamplingMethod
    # no params needed, the actual Household size will be inferred from the "Household" Setting in the Simulation.
end


#=
Now we have to define a second implementation of "sample_contacts()". Please note, that the first parameter has changed. The structs we defined here will be used to help the Julia Compiler to infer which implementation of "sample_contacts()" should be used.
To achieve "HouseholdSizeBasedContactSampling", we just need to fetch the current size of the individuals vector stored in the setting passed to this function. The parameter "setting" contains the setting of the individual, for which we perform the contact sampling.
=#
"""

# Parameters:
- `hh_size_based_sampling` = your own struct
- `setting` = setting of the `ego`
- `ego` = Individual, for which the contacts should be sampled
- `present_individuals` = Individuals, currently present in `setting` 
- `tick` = current tick of the simulation object
"""
function GEMS.sample_contacts(hh_size_based_sampling::HouseholdSizeBasedContactSampling, setting::Setting, ego::Individual, present_individuals::Vector{Individual}, tick::Int16)

    if isempty(present_individuals)
        throw(ArgumentError("No Individual is present in $setting. Please provide a Setting, where at least 1 Individual is present!"))
    end

    if length(present_individuals) == 1
        return Individual[]
    end

    num_of_contacts = size(setting) - 1
    res = Vector{Individual}(undef, num_of_contacts)
    cnt = 0
    # Draw until contact list is filled, skip whenever the index individual was selected
    while cnt < num_of_contacts
        contact = rand(present_individuals)
        # if contact is NOT index individual, add them to contact list
        if Ref(contact) .!== Ref(ego)
            res[cnt + 1] = contact
            cnt += 1
        end
    end

    return res
end


########### Case of batch runs lp=7 ip =5 alternative gs 17.5 alternative ###################################
# No  distributions for infectious period and onset ##################

rds = ResultData[]

for i in 1:100
    sim = nothing
    sim = Simulation(BASE_FOLDER * "/applications/workdays_weekends/weekends_paper/lp_7_ip_5_alt/synth_pop_lp7_ip5_alt.toml", label = "alternative latent period 7")

    run!(sim)
    rd = ResultData(sim, style = "LightRD")
    push!(rds, rd)
    println(i)
end

#rd = rds[1]
#gemsplot(rd, type = :TickCasesBySetting, xlims = (0, 500))

bd = BatchData(rds)

rep = buildreport(bd, "Report")
generate(rep, BASE_FOLDER  *  "/applications/workdays_weekends/weekends_paper/lp_7_ip_5_alt")


# get aggregated info on total infections
save_ag_data_base = tick_cases(bd)

jldsave(BASE_FOLDER*"/applications/workdays_weekends/weekends_paper/lp_7_ip_5_alt/infections_ag.jld2"; save_ag_data_base)
save_ag_data = JLD2.load(BASE_FOLDER  *"/applications/workdays_weekends/weekends_paper/lp_7_ip_5_alt/infections_ag.jld2")["save_ag_data_base"]

plot(save_ag_data.tick, [save_ag_data_base.mean], 
    #bar_width = 0.2,]
    xlabel = "day",
    ylabel = "Day cases",
    linecolor = :match,
    #xticks = data_wp_sizes,
label= ["mean baseline"])


########### END of batch runs ###################################


########### Case of 1 simulation ###################################
sim1 = nothing
sim2 = nothing
# onset 1
sim1 = Simulation(BASE_FOLDER * "/applications/workdays_weekends/synth_pop_05mln_wp100000_onset1_gs_inf_0.075.toml", label = "baseline")
sim2 = Simulation(BASE_FOLDER * "/applications/workdays_weekends/synth_pop_05mln_wp100000_onset1_gs_inf_0.075_alt.toml", label = "alternative")
run!(sim1)
run!(sim2)

rd1 = ResultData(sim1)
rd2 = ResultData(sim2)
gemsplot(rd1, type = :TickCasesBySetting, xlims = (0, 200))
gemsplot([rd1, rd2], type = :TickCasesBySetting, xlims = (0, 200))

gemsplot(rd1, type = :TickCases, xlims = (0, 200))
#gemsplot([rd1, rd2], type = :EffectiveReproduction, xlims = (0, 500))
#gemsplot([rd1 rd2], type = :LatencyHistogram)

gemsplot(rd1, type = :CumulativeCases, xlims = (0, 200))





comp = rd |> compartment_periods
mean(comp[!,"exposed"])
########### Case of 1 simulation END ###################################




################################ CALIBRATION CHECK OF SHARES ####################
infections_data = rd |> infections
total_inf = nrow(infections_data)
#vscodedisplay(infections_data)

#combine(groupby(infections_data, [:setting_type]), nrow => :count)
data_f = combine(groupby(infections_data, [:tick, :setting_type]), nrow => :count)

sort!(data_f)

#vscodedisplay(data_f)

res_settings_infections = combine(groupby(infections_data, [:setting_type]), nrow => :count)
sort!(res_settings_infections)
#vscodedisplay(res_settings_infections)
#-------------------------------------------------------------------
plt = bar(res_settings_infections.setting_type, res_settings_infections.count/total_inf, label= "Share of Infections by setting")

res_settings_infections.count/total_inf
################################ CALIBRATION CHECK OF SHARES ####################


onset = [1,2,3,4,5,6,7]
baseline = [8.6,	7.3, 6.9,	6.3,	6.2,	5.3,	4.9]
alternative = [5.4,	5.0,	6.1,	7.1,	8.1,	8.2,	5.1]


plt = plot(onset, [baseline, alternative], seriestype = :line, 
label = ["baseline" "alternative"], 
linestyle = :auto, 
marker = :circle,
xlabel = "symptoms onset period",
ylabel = "attack rate, %",
ylims = [-0.4,10])
