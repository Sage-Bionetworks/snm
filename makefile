cd ../../compiled/
rm -rf ./snm_1.0*
mkdir ./snm_1.0.0
cp -r ../SNM/trunk/ ./snm_1.0.0/
cd ./snm_1.0.0
rm README makefile
rm -rf .svn
cd ./man/
rm -rf .svn
rm buildBasisFunction.Rd edge.glr.Rd make.ref.model.matrices.Rd	sim.probe.specific.Rd buildBasisSplineMatrix.Rd	edge.qvalue.Rd make.snm.obj.Rd	calcArrayEffects.Rd err.msg.Rd makeDataObject.Rd fit.fast.model.Rd rm.zero.cols.Rd calculate.nulls.Rd fit.model.Rd snm.diagnosticPlot.Rd checkArguments.Rd get.coefs.Rd sim.function.var.Rd constructNormalizedData.Rd get.resids.Rd sim.intensity.dep.Rd edge.fit.Rd getSpanningSet.Rd
cd ../R
rm -rf .svn
rm checkArguments.R fit.fast.model.R get.coefs.R get.resids.R
cd ..
rm eme.R fix.lme4.R
rm -rf inst/.svn
rm -rf inst/doc/.svn
cd ..
read -p "Press any key to peform R CMD build, or press Control-C to end."
R CMD build snm_1.0.0
read -p "Press any key to peform R CMD check, or press Control-C to end."
R CMD check snm_1.0.0.tar.gz