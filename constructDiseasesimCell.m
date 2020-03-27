function [diseaseSim] = constructDiseasesimCell(meshSim,dtSim,dm, gamadd, gamall)
    DD1 = meshSim;
    DD2 = dtSim;
    [dmGauSim,~] = GaussianKernel(dm, gamadd, gamall);
    DD3 = dmGauSim;

    diseaseSim={DD1 DD2 DD3};

    clear DD1 DD2 DD3 dmGauSim dtSim meshSim

end

