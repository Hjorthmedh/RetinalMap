require(WhitelawCowan)

rgc.file=system.file("WhitelawCowan-2000-WT-retinal-gradients.txt",
  package='WhitelawCowan')
sc.file=system.file("WhitelawCowan-2000-WT-retinal-gradients.txt",
  package="WhitelawCowan")

op.file='WhitelawCowan-2000-WT-map.txt'

## We run for just a few iterations to check that everything can be run okay.

a <- run_WC(rgc.file,sc.file,op.file, stepbreaks=c(1,3,5))

#The output matrix is written to op.file in the current directory
