function [miRNASim] = constructMiRNAsimCell(seqSim,mtSim,meshSim,msSim,dm, gamadd, gamall)

    MM1 = seqSim;
    MM2 = mtSim;

    dmFunSim=dmFunMiRNASim(meshSim,dm);
    MM3=dmFunSim;

    [~,dmGauSim] = GaussianKernel(dm, gamadd, gamall);
    MM4 = dmGauSim;
    MM5 = msSim;
    miRNASim={MM1 MM2 MM3 MM4 MM5};
    clear MM1 MM2 MM3 MM4 MM5 seqSim mtSim dmFunSim dmGauSim msSim
end

