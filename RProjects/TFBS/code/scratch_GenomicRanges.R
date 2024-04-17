## -----------------------------------------------------------
## precede() and follow()
## -----------------------------------------------------------
query <- GRanges("A", IRanges(c(5, 20), width=1), strand="+")
#seqnames    ranges strand
#<Rle> <IRanges>  <Rle>
#  [1]        A         5      +
#  [2]        A        20      +
subject <- GRanges("A", IRanges(rep(c(10, 15), 2), width=1),
                   strand=c("+", "+", "-", "-"))
#seqnames    ranges strand
#<Rle> <IRanges>  <Rle>
#  [1]        A        10      +
#  [2]        A        15      +
#  [3]        A        10      -
#  [4]        A        15      -

#nämä riippuu juosteesta
precede(query, subject) #Mikään subject ei ole ennen queryä
#1 NA query[1] on ennen subjectia 1, query[2] ei ole ennen mitään
follow(query, subject)
# NA 2 query[2] seuraa jotakin subjectia[2], query[1] ei seuraa mitään
strand(query) <- "-"
precede(query, subject)
#NA  4 query[2] edeltää subjectia 4
follow(query, subject)
#3 NA

#Se on lähempänä joka on samassa strandissa
precede(query, subject, ignore.strand=TRUE) #Mikään subject ei ole ennen queryä
#1 NA query[1] on ennen subjectia 1, query[2] ei ole ennen mitään
follow(query, subject, ignore.strand=TRUE)
# NA 4 
strand(query) <- "-"
precede(query, subject, ignore.strand=TRUE)
#1 NA 
follow(query, subject, ignore.strand=TRUE)
#NA 4


## ties choose first in order
query <- GRanges("A", IRanges(10, width=1), c("+", "-", "*"))
subject <- GRanges("A", IRanges(c(5, 5, 5, 15, 15, 15), width=1),
                   rep(c("+", "-", "*"), 2))
precede(query, subject)
precede(query, rev(subject))

## ignore.strand=TRUE treats all ranges as '+'
precede(query[1], subject[4:6], select="all", ignore.strand=FALSE)
precede(query[1], subject[4:6], select="all", ignore.strand=TRUE)

## -----------------------------------------------------------
## nearest()
## -----------------------------------------------------------
## When multiple ranges overlap an "arbitrary" range is chosen
query <- GRanges("A", IRanges(5, 15))
subject <- GRanges("A", IRanges(c(1, 15), c(5, 19)))
nearest(query, subject)

## select="all" returns all hits
nearest(query, subject, select="all")

## Ranges in 'x' will self-select when 'subject' is present
query <- GRanges("A", IRanges(c(1, 10), width=5))
nearest(query, query)

## Ranges in 'x' will not self-select when 'subject' is missing
nearest(query)

## -----------------------------------------------------------
## nearestKNeighbors()
## -----------------------------------------------------------
## Without an argument, k defaults to 1
query <- GRanges("A", IRanges(c(2, 5), c(8, 15)))
subject <- GRanges("A", IRanges(c(1, 4, 10, 15), c(5, 7, 12, 19)))
nearestKNeighbors(query, subject)

## Return multiple neighbors with k > 1
nearestKNeighbors(query, subject, k=3)

## select="all" returns all hits
nearestKNeighbors(query, subject, select="all")

## -----------------------------------------------------------
## distance(), distanceToNearest()
## -----------------------------------------------------------
## Adjacent, overlap, separated by 1
query <- GRanges("A", IRanges(c(1, 2, 10), c(5, 8, 11)))
subject <- GRanges("A", IRanges(c(6, 5, 13), c(10, 10, 15)))
distance(query, subject)

## recycling
distance(query[1], subject)

## zero-width ranges
zw <- GRanges("A", IRanges(4,3))
stopifnot(distance(zw, GRanges("A", IRanges(3,4))) == 0L)
sapply(-3:3, function(i) 
  distance(shift(zw, i), GRanges("A", IRanges(4,3))))

query <- GRanges(c("chr1"), IRanges(c(10, 2), width=1),
strand=c("+", "-" ), name=c("TF1_HEAD_POS", "TF1_HEAD_NEG"))

subject <- GRanges("chr1", IRanges(c(5, 7), width=1), strand = c("+", "-"), name=c("TF2_TAIL_POS", "TF2_TAIL_NEG"))
#[1]        A         1      *
#[2]        B         5      *
distanceToNearest(query, subject)
distanceToNearest(query, subject, ignore.strand=TRUE)

distanceToNearest(rev(query), subject)
distanceToNearest(query, subject, ignore.strand=TRUE)


distanceToNearest(subject, query)
distanceToNearest(subject, query,ignore.strand=TRUE)


TF1_head=GRanges("chr1", IRanges(10), strand="-")
TF2_tail=GRanges("chr1", IRanges(10), strand="-")

distance(shift(TF1_head,1), TF2_tail) #0 just next to each other

distance(TF1_head, TF2_tail) #0 --> -1 one base pair overlap
start(TF1_head)-start(TF2_tail) #0

distance(shift(TF1_head,-1), TF2_tail) #0 --> -2 two base pair overlap
start(shift(TF1_head,-1))-start(TF2_tail) #-1

distance(shift(TF1_head,-2), TF2_tail) #1 --> -3 two base pair overlap
start(shift(TF1_head,-2))-start(TF2_tail) -2

distance(shift(TF1_head,-3), TF2_tail) #2 --> -4 two base pair overlap
start(shift(TF1_head,-3))-start(TF2_tail) -3


## distance() with GRanges and TxDb see the 
## ?'distance,GenomicRanges,TxDb-method' man 
## page in the GenomicFeatures package.