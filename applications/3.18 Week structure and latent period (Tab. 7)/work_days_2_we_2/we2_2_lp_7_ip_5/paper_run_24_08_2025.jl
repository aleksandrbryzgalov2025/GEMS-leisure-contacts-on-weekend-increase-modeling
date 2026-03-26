using GEMS, Parameters, DataFrames, TOML, Plots, FileIO, Distributions, CSV, CategoricalArrays, JLD2, Random
BASE_FOLDER = dirname(dirname(pathof(GEMS)))

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

@with_kw struct FixedContacts <: ContactSamplingMethod
    distribution::String = "Poisson"
    mean_number_of_contacts_weekday::Float64
    mean_number_of_contacts_weekend::Float64
    wdDistribution::Distribution = eval(Meta.parse(distribution))(mean_number_of_contacts_weekday...)
    weDistribution::Distribution = eval(Meta.parse(distribution))(mean_number_of_contacts_weekend...)
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
    if isempty(present_individuals)
        throw(ArgumentError("No Individual is present in $setting. Please provide a Setting, where at least 1 Individual is present!"))
    end

    if length(present_individuals) == 1
        return Individual[]
    end

    num_of_contacts = 0

    if (tick % 4) == 2 || (tick % 4) == 3  
        num_of_contacts = rand(fixed_contacts.weDistribution)
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



########### Case of batch runs lp= 7 ip=5, 2 after 2 WE days ###################################
#  ##################

rds = ResultData[]

for i in 1:100
    sim = nothing

# !!! check the folder path !!!

    sim = Simulation(BASE_FOLDER*"/applications/workdays_weekends/weekends_paper/we2_2_lp_7_ip_5/synth_pop_lp7_ip5_we2_2.toml", label = "Workdays 2 weekends 2 lp7")
    run!(sim)
    rd = ResultData(sim, style = "LightRD")
    push!(rds, rd)
    println(i)
end

bd = BatchData(rds)


rep = buildreport(bd, "Report")
# !!! check the folder path !!!
generate(rep, BASE_FOLDER  *  "/applications/workdays_weekends/weekends_paper/we2_2_lp_7_ip_5")


# get aggregated info on total infections
save_ag_data_base = tick_cases(bd)

# !!! check the folder path !!!
jldsave(BASE_FOLDER*"/applications/workdays_weekends/weekends_paper/we2_2_lp_7_ip_5/infections_ag.jld2"; save_ag_data_base)

# !!! check the folder path !!!
save_ag_data = JLD2.load(BASE_FOLDER  *"/applications/workdays_weekends/weekends_paper/we2_2_lp_7_ip_5/infections_ag.jld2")["save_ag_data_base"]

plot(save_ag_data.tick, [save_ag_data_base.mean], 
    #bar_width = 0.2,]
    xlabel = "day",
    ylabel = "Day cases",
    linecolor = :match,
    #xticks = data_wp_sizes,
label= ["mean baseline"])


########### END of batch runs ###################################


