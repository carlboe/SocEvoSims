# attempt to parametrize each simulation with an input file which will
# double as documentation.  This way, there is no way to introduct
# typos downstream, e.g. to use kinship5 to run the simulation but
# kinship3 in the analysis.

# the syntax of the definition file is that 


scan(file="simdefs.dat",what=list("","",""), comment.char="#",sep=":",

sim.num:Sim 141:  ## unique label for this simulation;     
                        ## if omitted, it will take the basename of
                        ## the working directory, e.g. sim141     
start.year:  1 : ## starting year of simulation
Ncycles:15000  : ## number of cycles to simulate
Ncheckpoint:500: ## checkpoints happen this often

mute.rate: 0.01:
additive.gene.risk: .1/5; :
