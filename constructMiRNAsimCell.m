function [miRNASim] = constructMiRNAsimCell(seqSim,mtSim,meshSim,dm, gamadd, gamall)

    MM1 = seqSim;
    MM2 = mtSim;

    dmFunSim=dmFunMiRNASim(meshSim,dm);
    MM3=dmFunSim;

    [~,dmGauSim] = GaussianKernel(dm, gamadd, gamall);
    MM4 = dmGauSim;

    miRNASim={MM1 MM2 MM3 MM4};
    clear MM1 MM2 MM3 MM4 seqSim mtSim dmFunSim dmGauSim
end

