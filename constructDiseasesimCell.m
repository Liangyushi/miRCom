function [diseaseSim] = constructDiseasesimCell(meshSim,dtSim,dsSim,dm,gamadd, gamall)
    DD1 = meshSim;
    DD2 = dtSim;
    [dmGauSim,~] = GaussianKernel(dm, gamadd, gamall);
    DD3 = dmGauSim;
    DD4=dsSim;

    diseaseSim={DD1 DD2 DD3 DD4};

    clear DD1 DD2 DD3 DD4 dmGauSim dtSim meshSim dsSim

end

