RetinalMap
==========

A framework for retinotopic map formation modelling and
evaluation. Currently four models are fully supported (Gierer, Koulakov,
Whitelaw and Willshaw). There is also code to run the Grimbert and
Cang model.

URL for github repository, this will create a RetinalMap directory:
https://github.com/Hjorthmedh/RetinalMap.git

##Installation

To be able to run all the models you need to compile the C
routines. They are located in the RetinalMap/@RetinalMap directory, and in the
base directory. The first one is for the Koulakov model, the second
one for the Willshaw model.
 
    mex stepFastGlobal.c
    mex stepMarkerInductionC.c 

The Gierer and Whitelaw models are distributed as R packages. The
detailed instructions for them can be found in a PDF in respective subdirectory.

In R:

    install.packages(c("sjevor"), contriburl="http://damtp.cam.ac.uk/user/eglen/r/")
    install.packages(c("fields"))

In the shell (make sure you are in the RetinalMap directory):

    R CMD INSTALL gierer
    R CMD INSTALL WhitelawCowan

Finally you need to install the Lattice analysis in the parent
directory of RetinalMap. To get the Lattice analysis code:

    git clone git://github.com/davidwillshaw/map-analysis.git

You should now have the two packages installed in the same directory, ie. YOUR_PATH/RetinalMap
and YOUR_PATH/map-analysis.

##Running the models

To run a simulation using the framework we use:

    report = ComparePhenotypeModelling(phenotype,action,expNum,plotFigures,kMask,model);

Here `phenotype` specifies which phenotype you want to simulate, your
standard options are `'WT'`,`'Isl2heterozygous'`, `'Isl2homozygous'`,
`'ephrinA2mm'`, `'ephrinA5mm'`, `'ephrinA2mmA5mm'` and `'TKO'`; `action` is
either `'run'` or `'analyse'`; `plotFigures` can be set to `true` or `false`;
`kMask` is normally `0`, it specifies if we should use any cis-interaction
to mask the gradients (not used in this study); `model` is one of `'Gierer2D'`, `'koulakov'`,
`'WhiteCow'` or `'Markerinduction'`.

For example:

    report = ComparePhenotypeModelling('WT','run',1,1,0,'koulakov');

Once the simulation finishes you can use the same command but with `'analyse'` as the action to get our first version of plots and summary sheets. However there is another version which has replaced it:

    report = makeFiguresForCompareArticle(phenotype,expRange,model,plotFigures)

In our earlier example:

    report = makeFiguresForCompareArticle('WT',1,'koulakov',1)

I would recommend using this second option, as it generates the figures relevant to the paper.


##Generating the figures in the paper

If you want to regenerate all the modelling figures for the paper there are two commands you need to do.

To run all simulations (called from matlab)

    runCompareArticleKoulakov
    runCompareArticleGierer
    runCompareArticleWhitelaw
    runCompareArticleWillshaw

To analyse all simulations:
    runCompareArticleMacAnalyseAll

This will output the figures to `FIGS/ComparePhenotype`


##Running your first own simulation

If you want more freedome in defining your own simulations, then you
can work directly with the `RetinalMap` class.

    r  = RetinalMap('experiments/myfirstexperiment.txt');
    r.run();
	r.plotMapForward(1);

For more instructions see the `readme.rtf` file.
