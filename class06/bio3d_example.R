require(bio3d)

prot_interact <- function(pdb.name, chain = "A", elety = "CA", ret = F) {
  pdb.dat <- read.pdb(pdb.name)
  
  # select the particular chain and only the element type
  chain <- trim.pdb(pdb.dat, chain = chain, elety = elety)
  chain.b <- chain$atom$b
  
  plotb3(chain.b, sse = chain, typ = 'l', ylab = 'Bfactor')
  
  if(ret) {
    return(chain.b)
  }
}