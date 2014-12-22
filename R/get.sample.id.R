
.get.sample.id <- function(PARAMETERS) {
# number of samples
no_ip=length(PARAMETERS$IP_BAM)
no_input=length(PARAMETERS$INPUT_BAM)
no_treated_ip=length(PARAMETERS$TREATED_IP_BAM)
no_treated_input=length(PARAMETERS$TREATED_INPUT_BAM)

# id
id_untreated_ip=1:no_ip
id_untreated_input=(no_ip+1):(no_ip+no_input)
if (no_treated_ip >0) {
id_treated_ip=(no_ip+no_input+1):(no_ip+no_input+no_treated_ip) } else {id_treated_ip = numeric(0)}

if (no_treated_input > 0) {
id_treated_input=(no_ip+no_input+no_treated_ip+1):(no_ip+no_input+no_treated_ip+no_treated_input) } else {id_treated_input=numeric(0)}

# together id
id_ip=c(id_untreated_ip,id_treated_ip)
id_input=c(id_untreated_input,id_treated_input)

# sample names
sample_names=c(rep("untreated_ip",no_ip),rep("untreated_input",no_input),
               rep("treated_ip",no_treated_ip),rep("treated_input",no_treated_input))

# result
id=list(ip=id_ip,input=id_input,untreated_ip=id_untreated_ip,untreated_input=id_untreated_input,
        treated_ip=id_treated_ip,treated_input=id_treated_input,sample_names=sample_names)

return(id)
}