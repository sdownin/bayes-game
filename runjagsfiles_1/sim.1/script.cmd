model in "model.txt"
data in "data.txt"
compile, nchains(4)
parameters in "inits1.txt", chain(1)
parameters in "inits2.txt", chain(2)
parameters in "inits3.txt", chain(3)
parameters in "inits4.txt", chain(4)
initialize
adapt 267
update 1067
monitor q, thin(1)
monitor epsilon, thin(1)
update 2667
parameters to "out1.Rdump", chain(1)
parameters to "out2.Rdump", chain(2)
parameters to "out3.Rdump", chain(3)
parameters to "out4.Rdump", chain(4)
coda *, stem(sim.1/CODA)
samplers to sim.1/samplers.csv
update 0
model clear
exit
