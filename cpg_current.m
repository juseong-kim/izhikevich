function I = cpg_current(tstop,dt,start,iLM,iRM,iLC,iRC,iLL,iRL,iLE,iRE,iLR,iRR)
    si = start / dt;
    I = [zeros(1,si-1), iLM * ones(1,(tstop + 1)/dt-si);
            zeros(1,si-1), iRM * ones(1,(tstop + 1)/dt-si);

            zeros(1,si-1), iLC * ones(1,(tstop + 1)/dt-si);
            zeros(1,si-1), iRC * ones(1,(tstop + 1)/dt-si);

            zeros(1,si-1), iLL * ones(1,(tstop + 1)/dt-si);
            zeros(1,si-1), iRL * ones(1,(tstop + 1)/dt-si);

            zeros(1,si-1), iLE * ones(1,(tstop + 1)/dt-si);
            zeros(1,si-1), iRE * ones(1,(tstop + 1)/dt-si);

            zeros(1,si-1), iLR * ones(1,(tstop + 1)/dt-si);
            zeros(1,si-1), iRR * ones(1,(tstop + 1)/dt-si)
        ];
end