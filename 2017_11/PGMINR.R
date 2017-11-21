# Code from 
# http://people.math.aau.dk/~sorenh/misc/2012-Oslo-GMwR/GMwR-notes.pdf


require(gRbase)
require(gRain)
require(gRim)
require(Rgraphviz)
require(RBGL)

########################################################
# LET'S LOOK AT SYNTAX FOR CONSTRUCTING GRAPHS
########################################################

# UNDIRECTED GRAPHS
# use function "ug" to create undirected graph
g1 <- ug(~a:b:e + a:c:e + b:d:e + c:d:e + c:g + d:f)
class(g1)
as(g1, "matrix") # notice matrix is symmetric because undirected
plot(g1)

# can also create an undirected graph with matrices
m <- as(g1, "matrix") 
as(m, "graphNEL")
plot(as(m, "graphNEL"))

g1a <- addEdge("a", "d", g1)
g1b <- removeEdge("c", "d", g1)
par(mfrow = c(1, 3))
plot(g1, main = "g1")
plot(g1a, main = "g1a")
plot(g1b, main = "g1b")

g1c <-subGraph(c("b", "c", "d", "e"), g1)
par(mfrow = c(1, 2))
plot(g1, main = "g1")
plot(g1c, main = "g1c")


########################################################
# LET'S LOOK AT OPERATIONS FOR UNDIRECTED GRAPHS
########################################################

## a clique is a maximal complete subset not contained in a larger complete subset
plot(g1)
is.complete(g1)
is.complete(g1, set = c("a", "b", "e"))
is.complete(g1, set = c("a", "b"))
is.complete(g1, set = c("f", "g"))
is.complete(g1, set = c("f"))

# but we don't want to try this all out by hand...
str(maxClique(g1))


# a subset S is said to separate subsets A and B if every path between a vertex in A and a vertex in B 
# contains a vertex from S
g2 <- ug(~a:b:e + a:c:e + b:d:e + c:d:e)
plot(g2)
separates("a", "d", c("b", "c", "e"), g2)
separates("a", "d", c("b", "e"), g2)
separates(c("a" ,"b"), "d", c("e"), g2)

# If A and B are separated by S in the dependency graph then A and B are conditionally independent given S
plot((g3 <- ug(~ A:B + B:C:D + C:E + D:E)))
separates(c("D", "E"), "A", "B", g3)
separates(c("D", "E"), "A", c("B", "C"), g3)

########################################################
# LET'S LOOK AT CONSTRUCTING DIRECTED GRAPHS
########################################################
# DIRECTED GRAPHS
# create from dag function in gRbase
# can specify by a list of formulas or a list of vectors
dag0 <- dag(~a, ~b * a, ~c * a * b, ~d * c * e, ~e * a)
plot(dag0)
dag0 <- dag(~a + b:a + c:a:b + d:c:e + e:a)
plot(dag0)
dag0 <- dag("a", c("b", "a"), c("c", "a", "b"), c("d", "c", "e"), c("e", "a"))
plot(dag0)
dag0
# can also create from a matrix, but now there is directionality AND MATRX SHOULD NOT BE SYMMETRIC

########################################################
# OPERATIONS WITH DIRECTED GRAPHS
########################################################

## MORALIZATION
## to moralize a graph (1) add edges between the parents and each node and then (2) replace all directed edges with undirected ones, thus returning an undirected graph
dagM <- moralize(dag0)
# why is this called moralization?

# INDEPENDENCES
## ANCESTRAL GRAPHS
## if there is a path from a to b, we write a->b The ancestors of a node b are the nodes a such that a-> b. The ancestral set of a set A is the union of A with its ancestors. The ancestral graph of a set A is the subgraph induced by the ancestral set of A
parents("d", dag0)
children("c", dag0)

par(mfrow = c(1, 2))
ancestralSet(c("a", "c", "e"), dag0)
plot(ancestralGraph(c("a", "c", "e"), dag0))
plot(dag0)

# to check if A and B are conditionally independent on S from the ancestral graph of A union B union S. Moralize this ancestral graph. If A and B are separated by S in this moral graph, then A indp B given S
plot(moralize(ancestralGraph(c("a", "c", "e"), dag0)))
testG <- moralize(ancestralGraph(c("a", "c", "e"), dag0))
separates("c", "e", "a", testG)

########################################################
# COMPUTING MARGINAL DISTRIBUTIONS
########################################################

# SMALL WORKED EXAMPLE
# F = FEVER
# T = TEMPERATURE
# H = HEADACHE
plot((FTH <- dag(~ F + T:F + H:T)), "circo")

(p.F <- parray("F", levels = 2, values = c(0.01, 0.99)))

# t given f
(p.TgF <- parray(c("T", "F"), levels = c(2, 2), values = c(0.95, 0.05, 0.001,0.999)))

# h given t
(p.HgT <- parray(c("H", "T"), levels = c(2, 2), values = c(0.80, 0.20, 0.01, 0.99)))

# brute force computations
# calculate joint distribution
p.FT <- tableMult(p.F, p.TgF)
p.FTH <- tableMult(p.FT, p.HgT)
as.data.frame.table(p.FTH)

# calculate marginal distribution p FH (page 20 of tutorial)
p.FH <- tableMargin(p.FTH, margin=c("F", "H"))
as.data.frame.table(p.FH)

# calculate conditional distribution p(F|H)
p.H <- tableMargin(p.FH, margin = 'H')
(p.FgH <- tableDiv(p.FH, p.H))

# always want to avoid calculating the joint distribution
plot(FTH)
plot(moralize(FTH))

# clique potential representation p(F,T,H) = f(F,T) * f(T, H) / f(T)
(qFT <- tableMult(p.F, p.TgF))
(qTH <- p.HgT)
(qT <- parray("T", levels = 2, values = 1))

# working inwards in junction tree
# work inwards towards root, from FT towards TH
# fT* = sum over F of f(F,T)
# f(T,H)* = f(T,H) * f(T)*/f(T)
(qTs <- tableMargin(qFT, "T"))
(qTHs <- tableMult(qTH, tableDiv(qTs, qT)))

# now iterate back working 'outward' in junction tree
# working outwards in junction tree
# work outwards from root (from TH to FT)
(qTss <- tableMargin(qTHs, "T"))
(qFTs <- tableMult(qFT, tableDiv(qTss, qTs)))

# this leaves us with marginal distributions on all cliques and separators
qFTs
qTHs
qTs

# probability of temperature
qTs

# probability of fever
tableMargin(qFT, "F")

# probability of headache
tableMargin(qTH, "H")

# propagating findings
# support we have the finding H = yes = H1
# set any entry in f(T,H) inconsistent with H = H1 to 0 and this yields a new potential
# repeat computations above
qTH
## set the finding H = H1
qTH[c(2,4)] <- 0
qTH

# repeat everything that was done above
(qTs <- tableMargin(qFT, "T"))
(qTHs <- tableMult(qTH, tableDiv(qTs, qT))) # <----using the new qTH
(qTss <- tableMargin(qTHs, "T"))
(qFTs <- tableMult(qFT, tableDiv(qTss, qTs)))


# but the table only contains the clique probabilities up to a normalizing constant
sum(qFTs)

# to get probability of a fever we must normalize
tableMargin(qFTs, "F")/sum(qFTs)

# after working inwards and outwards in the junction tree, clique potentials are consistent - they match
# on their separators
tableMargin(qFTs, "T")
tableMargin(qTHs, "T")
qTss

########################################################
# REALISTIC EXAMPLE
########################################################

# OF COURSE THIS ALL COMES MORE WRAPPED UP AS WELL
# 1 minute guide to gRain
# https://cran.r-project.org/web/packages/gRain/vignettes/gRain-intro.pdf
# now let's look at higher level functionality
yn <- c("yes","no")
a <- cptable(~asia, values=c(1,99),levels=yn) 
t.a <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn) 
s <- cptable(~smoke, values=c(5,5), levels=yn) 
l.s <- cptable(~lung|smoke, values=c(1,9,1,99), levels=yn) 
b.s <- cptable(~bronc|smoke, values=c(6,4,3,7), levels=yn) 
e.lt <- cptable(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn) 
x.e <- cptable(~xray|either, values=c(98,2,5,95), levels=yn) 
d.be <- cptable(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)

plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
plist
plist$either ## Notice: a logical node
net1 <- grain(plist)
plot(net1)

# the network can be queried to give marginal probabilities
querygrain(net1, nodes=c("lung","bronc"), type="marginal")

# can also query for a joint distribution (brute force)
querygrain(net1,nodes=c("lung","bronc"), type="joint")

# can enter hard evidence
net12 <- setEvidence(net1, evidence=list(asia="yes", dysp="yes")) 
net12 <- setEvidence(net1, + nodes=c("asia", "dysp"), states=c("yes", "yes"))

# probability of observing the evidence
pEvidence(net12)

# can query the network given the evidence
querygrain( net12, nodes=c("lung","bronc") )
querygrain( net12, nodes=c("lung","bronc"), type="joint" )

# specifying hard evidence
# notice we can reuse cpt specified above
(plist1 <- compileCPT(list(a, t.a)))
(chest1 <- grain(plist1))
querygrain(chest1)

setEvidence(chest1, evidence = list(asia = "yes"))
querygrain(setEvidence(chest1, evidence = list(asia = "yes")))

# specifying virtual evidence (likelihood evidence)
# 2 ways to do this
# first, more explicit way is to build a network with a new node
g.a <- parray(c("guess.asia", "asia"), levels = list(yn, yn), values = c(0.8, 0.2, 0.1, 0.9))
plist2 <- compileCPT(list(a, ta., g.a))
chest2 <- grain(plist2)
querygrain(chest2)
querygrain(setEvidence(chest2, evidence=list(guess.asia = "yes")))

# second way is can pass in likelihood evidence directly for original graph
querygrain(setEvidence(chest1, evidence = list(asia = c(1,0))))











